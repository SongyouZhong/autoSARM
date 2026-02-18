"""
SARM 异步任务处理器
支持并发处理和进度更新
使用 PostgreSQL 数据库和 SeaweedFS 对象存储
"""

import asyncio
import json
import logging
import os
import shutil
import time
import uuid
from argparse import Namespace
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Optional, Dict, Any

import asyncpg

from autosarm.config.settings import settings, get_default_cores
from autosarm.storage import get_storage

logger = logging.getLogger("async_task_processor")


class TaskProgressCallback:
    """任务进度回调类"""

    def __init__(self, task_id: str, connection, processor):
        self.task_id = task_id
        self.connection = connection
        self.processor = processor
        self._is_completed = False

    async def update_progress(self, progress: float, info: str = None):
        """更新任务进度"""
        if self._is_completed:
            return

        try:
            self.processor.task_progress[self.task_id] = {
                "overall_progress": progress,
                "current_step": info or "Processing...",
                "status": "processing",
                "last_updated": time.time(),
            }

            await self.connection.execute(
                "UPDATE tasks SET status = $1 WHERE id = $2",
                "processing", self.task_id
            )

            logger.info("Task %s progress: %.1f%% - %s",
                         self.task_id, progress, info or "")
        except Exception as e:
            logger.error("Failed to update progress for task %s: %s", self.task_id, e)

    def mark_completed(self):
        """标记任务为已完成"""
        self._is_completed = True


class AsyncTaskProcessor:
    """
    异步任务处理器
    
    支持多容器 Worker 模式:
    - 使用数据库行级锁 (SELECT FOR UPDATE SKIP LOCKED) 防止任务重复处理
    - 每个实例每次只处理一个任务 (max_workers=1)
    - 支持 docker compose 启动多个实例水平扩展
    """

    def __init__(self):
        task_settings = settings().task_processor
        db_settings = settings().database

        self.poll_interval = task_settings.poll_interval
        self.active_tasks: Dict[str, asyncio.Task] = {}
        self.task_progress: Dict[str, Dict[str, Any]] = {}
        self.is_running = True
        self.polling_task = None
        self._db_pool: Optional[asyncpg.Pool] = None

        # 每个容器实例每次只处理一个任务，便于水平扩展
        self.max_workers = 1
        self.thread_executor = ThreadPoolExecutor(max_workers=self.max_workers)

        # 数据库配置
        self.db_config = {
            'host': db_settings.host,
            'port': db_settings.port,
            'user': db_settings.user,
            'password': db_settings.password,
            'database': db_settings.database,
        }

        # 生成唯一的 worker ID 用于日志追踪
        self.worker_id = str(uuid.uuid4())[:8]

        logger.info("AsyncTaskProcessor initialized (worker_id=%s, max_workers=%d)",
                     self.worker_id, self.max_workers)

    async def start_polling(self):
        """启动数据库轮询"""
        if self.polling_task is None:
            await self._init_db_pool()
            logger.info("Starting database polling for sarm_analysis tasks...")
            self.polling_task = asyncio.create_task(self._poll_database_tasks())

    async def _init_db_pool(self):
        """初始化数据库连接池"""
        if self._db_pool is None:
            try:
                logger.info("Creating PostgreSQL connection pool...")
                self._db_pool = await asyncpg.create_pool(
                    host=self.db_config['host'],
                    port=self.db_config['port'],
                    user=self.db_config['user'],
                    password=self.db_config['password'],
                    database=self.db_config['database'],
                    min_size=1,
                    max_size=5,
                )
                logger.info("PostgreSQL connection pool created successfully")
            except Exception as e:
                logger.error(f"Failed to create database pool: {e}")
                raise

    async def _poll_database_tasks(self):
        """
        定时从数据库获取待处理的 SARM 分析任务
        
        使用 SELECT FOR UPDATE SKIP LOCKED 实现行级锁:
        - 防止多个 worker 同时获取同一个任务
        - SKIP LOCKED 确保如果任务被锁定，则跳过而不是等待
        - 每个 worker 每次只获取一个任务 (单任务模式)
        """
        while self.is_running:
            try:
                # 如果已有任务在执行，则等待
                if len(self.active_tasks) >= self.max_workers:
                    logger.debug("[Worker %s] Already processing %d task(s), waiting...",
                                 self.worker_id, len(self.active_tasks))
                    await asyncio.sleep(self.poll_interval)
                    continue

                async with self._db_pool.acquire() as connection:
                    # 使用事务和行级锁获取任务
                    async with connection.transaction():
                        task = await connection.fetchrow(
                            """
                            SELECT id, job_dir 
                            FROM tasks 
                            WHERE task_type = 'sarm_analysis' 
                              AND status = 'pending' 
                            ORDER BY created_at ASC
                            LIMIT 1
                            FOR UPDATE SKIP LOCKED
                            """
                        )

                        if task:
                            task_id, job_dir = task['id'], task['job_dir']

                            if task_id not in self.active_tasks:
                                # 立即将任务状态更新为 processing，防止其他 worker 获取
                                await connection.execute(
                                    "UPDATE tasks SET status = $1, started_at = NOW() WHERE id = $2",
                                    "processing", task_id
                                )

                                logger.info("[Worker %s] Claimed task: %s",
                                            self.worker_id, task_id)

                                # 提交任务到执行队列
                                await self.submit_task(task_id, job_dir)
                            else:
                                logger.debug("[Worker %s] Task %s already in progress, skipping",
                                             self.worker_id, task_id)
                        else:
                            logger.debug("[Worker %s] No pending tasks available", self.worker_id)

            except Exception as e:
                logger.error("[Worker %s] Error polling database for tasks: %s",
                             self.worker_id, e)

            await asyncio.sleep(self.poll_interval)

        logger.info("[Worker %s] Database polling stopped", self.worker_id)

    async def get_db_connection(self):
        """从连接池获取数据库连接"""
        try:
            if self._db_pool is None:
                await self._init_db_pool()
            return await self._db_pool.acquire()
        except Exception as e:
            logger.error(f"数据库连接失败: {e}")
            return None

    async def release_db_connection(self, connection):
        """释放数据库连接回连接池"""
        if connection and self._db_pool:
            await self._db_pool.release(connection)

    def _get_temp_dir(self) -> Path:
        """获取临时目录"""
        storage_settings = settings().storage
        temp_dir = Path(storage_settings.temp_dir)
        temp_dir.mkdir(parents=True, exist_ok=True)
        return temp_dir

    def _is_seaweedfs_path(self, job_dir: str) -> bool:
        """
        判断 job_dir 是否是 SeaweedFS 存储前缀格式
        
        本地路径格式: /tmp/astramolecula/jobs/sarm_analysis/{job_id}
        SeaweedFS 前缀格式: jobs/sarm_analysis/{job_id}
        """
        if job_dir.startswith('/'):
            return False
        return True

    def _convert_to_storage_prefix(self, job_dir: str) -> str:
        """
        将 job_dir 转换为 SeaweedFS 存储前缀
        
        本地路径: /tmp/astramolecula/jobs/sarm_analysis/{job_id}
        转换为:   jobs/sarm_analysis/{job_id}
        """
        if self._is_seaweedfs_path(job_dir):
            return job_dir

        parts = job_dir.split('/tmp/astramolecula/')
        if len(parts) > 1:
            return parts[1]

        if '/jobs/' in job_dir:
            idx = job_dir.index('/jobs/') + 1
            return job_dir[idx:]

        logger.warning("Cannot convert job_dir to storage prefix: %s", job_dir)
        return job_dir

    async def _fetch_task_params(self, connection, task_id: str) -> Dict[str, Any]:
        """
        从 sarm_task_params 表获取任务参数
        
        Args:
            connection: 数据库连接
            task_id: 任务ID
            
        Returns:
            任务参数字典
        """
        row = await connection.fetchrow(
            """
            SELECT task_subtype, csv_filename, analysis_type, value_columns,
                   log_transform, minimum_site1, minimum_site2, n_jobs, csv2excel,
                   fragment_core, root_title, input_file, tree_content,
                   highlight_dict, max_level
            FROM sarm_task_params
            WHERE task_id = $1
            """,
            task_id
        )

        if not row:
            logger.warning("No sarm_task_params found for task %s, using defaults", task_id)
            return {"task_subtype": "sarm"}

        params = dict(row)

        # 解析 JSON 字段
        for json_field in ('value_columns', 'tree_content', 'highlight_dict'):
            if params.get(json_field):
                try:
                    params[json_field] = json.loads(params[json_field])
                except (json.JSONDecodeError, TypeError):
                    pass

        logger.info("Loaded task params for %s: subtype=%s", task_id, params.get('task_subtype'))
        return params

    async def _download_input_files(self, storage_prefix: str, temp_input_dir: Path,
                                     csv_filename: str = "compounds.csv") -> None:
        """
        从 SeaweedFS 下载输入文件到临时目录
        
        Args:
            storage_prefix: SeaweedFS 存储前缀
            temp_input_dir: 本地临时输入目录
            csv_filename: CSV 文件名
        """
        storage = get_storage()

        # 下载 CSV 输入文件
        csv_key = f"{storage_prefix}/input/{csv_filename}"
        csv_local = temp_input_dir / csv_filename
        try:
            await storage.download_file(csv_key, csv_local)
            logger.info("Downloaded CSV file: %s", csv_key)
        except Exception as e:
            logger.error("Failed to download CSV file: %s", e)
            raise FileNotFoundError(f"CSV file not found in storage: {csv_key}")

    async def _upload_results_to_storage(self, task_id: str, job_dir: str,
                                          storage_prefix: str = None):
        """
        上传任务结果到 SeaweedFS
        
        Args:
            task_id: 任务ID
            job_dir: 本地任务目录（临时目录）
            storage_prefix: SeaweedFS 存储前缀
        """
        try:
            storage = get_storage()
            output_dir = Path(job_dir) / "output"

            if not output_dir.exists():
                logger.warning("Output directory not found: %s", output_dir)
                return

            uploaded_count = 0
            for file_path in output_dir.rglob("*"):
                if file_path.is_file():
                    relative_path = file_path.relative_to(output_dir)

                    if storage_prefix:
                        remote_key = f"{storage_prefix}/output/{relative_path}"
                    else:
                        remote_key = f"tasks/{task_id}/sarm/output/{relative_path}"

                    try:
                        await storage.upload_file(file_path, remote_key)
                        uploaded_count += 1
                    except Exception as e:
                        logger.error("Failed to upload file %s: %s", file_path, e)

            logger.info("Uploaded %d files to SeaweedFS for task %s", uploaded_count, task_id)

        except Exception as e:
            logger.error("Failed to upload results to storage for task %s: %s", task_id, e)

    def _run_sarm_sync(self, args: Namespace) -> None:
        """同步执行 SARM 矩阵生成（在线程池中运行）"""
        from autosarm.cli.create_sarm import run_sarm
        run_sarm(args)

    def _run_tree_sync(self, args: Namespace) -> None:
        """同步执行 SAR 树生成（在线程池中运行）"""
        from autosarm.cli.create_tree import run_tree
        run_tree(args)

    async def process_sarm_task(self, task_id: str, job_dir: str):
        """
        处理 SARM 分析任务（支持 SeaweedFS 存储）
        
        流程:
        1. 从数据库获取任务参数 (sarm_task_params)
        2. 从 SeaweedFS 下载输入文件
        3. 根据 task_subtype 执行 sarm 或 tree 计算
        4. 上传结果到 SeaweedFS
        5. 更新任务状态
        """
        connection = None
        temp_job_dir = None

        try:
            connection = await self.get_db_connection()
            if not connection:
                raise Exception("Failed to connect to database")

            progress_callback = TaskProgressCallback(task_id, connection, self)
            await progress_callback.update_progress(0, "Starting SARM analysis")

            original_cwd = os.getcwd()

            try:
                # 获取任务参数
                params = await self._fetch_task_params(connection, task_id)
                task_subtype = params.get('task_subtype', 'sarm')

                # 将 job_dir 转换为 SeaweedFS 存储前缀
                storage_prefix = self._convert_to_storage_prefix(job_dir)
                logger.info("Task %s: storage_prefix=%s (original job_dir=%s)",
                            task_id, storage_prefix, job_dir)

                # 创建临时目录
                temp_base = self._get_temp_dir()
                temp_job_dir = temp_base / task_id
                temp_job_dir.mkdir(parents=True, exist_ok=True)
                temp_input_dir = temp_job_dir / "input"
                temp_input_dir.mkdir(exist_ok=True)
                temp_output_dir = temp_job_dir / "output"
                temp_output_dir.mkdir(exist_ok=True)

                logger.info("Task %s: Created temp directory: %s", task_id, temp_job_dir)

                # 根据子类型分发执行
                if task_subtype == 'sarm':
                    await self._process_sarm_subtype(
                        task_id, params, storage_prefix,
                        temp_input_dir, temp_output_dir, temp_job_dir,
                        progress_callback
                    )
                elif task_subtype == 'tree':
                    await self._process_tree_subtype(
                        task_id, params, storage_prefix,
                        temp_input_dir, temp_output_dir, temp_job_dir,
                        progress_callback
                    )
                else:
                    raise ValueError(f"Unknown task_subtype: {task_subtype}")

                # 上传结果
                await progress_callback.update_progress(90, "Uploading results to storage")
                await self._upload_results_to_storage(task_id, str(temp_job_dir), storage_prefix)

                # 更新任务状态为完成
                await connection.execute(
                    "UPDATE tasks SET status = $1, finished_at = NOW() WHERE id = $2",
                    "finished", task_id
                )

                progress_callback.mark_completed()
                logger.info("Task %s completed successfully", task_id)

                self.task_progress[task_id] = {
                    "overall_progress": 100,
                    "current_step": "Completed",
                    "details": "SARM analysis completed successfully",
                    "status": "finished",
                    "last_updated": time.time(),
                }

            finally:
                os.chdir(original_cwd)

                # 清理临时目录
                if temp_job_dir and temp_job_dir.exists():
                    try:
                        shutil.rmtree(temp_job_dir, ignore_errors=True)
                        logger.info("Task %s: Cleaned up temp directory: %s", task_id, temp_job_dir)
                    except Exception as e:
                        logger.warning("Task %s: Failed to cleanup temp directory: %s", task_id, e)

        except Exception as e:
            logger.error("Task %s failed: %s", task_id, str(e))

            self.task_progress[task_id] = {
                "overall_progress": 0,
                "current_step": "Failed",
                "details": f"Task failed: {str(e)}",
                "status": "failed",
                "last_updated": time.time(),
            }

            if connection:
                try:
                    await connection.execute(
                        "UPDATE tasks SET status = $1, finished_at = NOW() WHERE id = $2",
                        "failed", task_id
                    )
                except Exception as db_error:
                    logger.error("Failed to update task status in database: %s", db_error)

        finally:
            if connection:
                await self.release_db_connection(connection)

            if task_id in self.active_tasks:
                del self.active_tasks[task_id]

    async def _process_sarm_subtype(self, task_id: str, params: Dict[str, Any],
                                     storage_prefix: str,
                                     temp_input_dir: Path, temp_output_dir: Path,
                                     temp_job_dir: Path,
                                     progress_callback: TaskProgressCallback):
        """
        处理 SARM 矩阵生成子任务
        
        Args:
            task_id: 任务ID
            params: 任务参数
            storage_prefix: SeaweedFS 存储前缀
            temp_input_dir: 临时输入目录
            temp_output_dir: 临时输出目录
            temp_job_dir: 临时任务根目录
            progress_callback: 进度回调
        """
        csv_filename = params.get('csv_filename', 'compounds.csv')

        # 下载输入文件
        await progress_callback.update_progress(5, "Downloading input files from storage")
        await self._download_input_files(storage_prefix, temp_input_dir, csv_filename)

        csv_file = temp_input_dir / csv_filename
        if not csv_file.exists():
            raise FileNotFoundError(f"CSV file not found: {csv_file}")

        await progress_callback.update_progress(10, "Validating input files")

        # 解析 value_columns
        value_columns = params.get('value_columns', [])
        if isinstance(value_columns, str):
            try:
                value_columns = json.loads(value_columns)
            except (json.JSONDecodeError, TypeError):
                value_columns = [value_columns]
        if not value_columns:
            # 尝试自动检测活性列
            import pandas as pd
            df_preview = pd.read_csv(csv_file, nrows=5)
            # 排除常见的非活性列
            exclude_cols = {'smiles', 'SMILES', 'Cano_SMILES', 'Name', 'name', 'ID', 'id',
                            'Scaffold', 'scaffold', 'Compound', 'compound'}
            numeric_cols = [c for c in df_preview.select_dtypes(include=['number']).columns
                           if c not in exclude_cols]
            if numeric_cols:
                value_columns = numeric_cols
                logger.info("Auto-detected value columns: %s", value_columns)
            else:
                raise ValueError("No value columns specified and auto-detection failed")

        # 构造 argparse Namespace 对象
        save_folder = str(temp_output_dir / "SAR_Results")
        # 安全获取参数，处理 None 值（数据库 NULL）
        n_jobs_raw = params.get('n_jobs')
        n_jobs = int(n_jobs_raw) if n_jobs_raw is not None else get_default_cores()

        minimum_site1_raw = params.get('minimum_site1')
        minimum_site1 = float(minimum_site1_raw) if minimum_site1_raw is not None else 3.0

        minimum_site2_raw = params.get('minimum_site2')
        minimum_site2 = float(minimum_site2_raw) if minimum_site2_raw is not None else 3.0

        args = Namespace(
            csvFile=str(csv_file),
            type=params.get('analysis_type', 'smiles') or 'smiles',
            column=value_columns,
            log=1 if params.get('log_transform', False) else 0,
            minimumSite1=minimum_site1,
            minimumSite2=minimum_site2,
            n_jobs=n_jobs,
            save_folder=save_folder,
            csv2excel=1 if params.get('csv2excel', False) else 0,
        )

        await progress_callback.update_progress(20, "Running SARM matrix generation")
        logger.info("Task %s: Starting SAR matrix generation with args: %s", task_id, vars(args))

        # 在线程池中运行 CPU 密集型计算
        loop = asyncio.get_event_loop()
        await loop.run_in_executor(self.thread_executor, self._run_sarm_sync, args)

        await progress_callback.update_progress(85, "SARM matrix generation complete")

    async def _process_tree_subtype(self, task_id: str, params: Dict[str, Any],
                                     storage_prefix: str,
                                     temp_input_dir: Path, temp_output_dir: Path,
                                     temp_job_dir: Path,
                                     progress_callback: TaskProgressCallback):
        """
        处理 SAR 树生成子任务
        
        Args:
            task_id: 任务ID
            params: 任务参数
            storage_prefix: SeaweedFS 存储前缀
            temp_input_dir: 临时输入目录
            temp_output_dir: 临时输出目录
            temp_job_dir: 临时任务根目录
            progress_callback: 进度回调
        """
        # SAR 树需要先前的 SARM 矩阵结果作为输入
        # 从 SeaweedFS 下载已有的 SAR_Results 目录
        await progress_callback.update_progress(5, "Downloading SAR results from storage")

        storage = get_storage()
        input_file = params.get('input_file', 'input.csv')

        # 下载 input.csv（树生成的输入数据）
        csv_key = f"{storage_prefix}/input/{input_file}"
        csv_local = temp_input_dir / input_file
        try:
            await storage.download_file(csv_key, csv_local)
        except FileNotFoundError:
            logger.warning("Input file not found: %s, trying alternative path", csv_key)

        # 下载先前 SARM 生成的结果（如果有）
        sar_results_key = f"{storage_prefix}/output/SAR_Results"
        sar_results_dir = temp_output_dir / "SAR_Results"
        sar_results_dir.mkdir(parents=True, exist_ok=True)
        try:
            await storage.download_directory(sar_results_key, sar_results_dir)
            logger.info("Downloaded existing SAR results")
        except Exception:
            logger.info("No existing SAR results to download")

        # 解析 tree_content 和 highlight_dict
        tree_content = params.get('tree_content', ['double-cut'])
        if isinstance(tree_content, str):
            try:
                tree_content = json.loads(tree_content)
            except (json.JSONDecodeError, TypeError):
                tree_content = ['double-cut']

        highlight_dict = params.get('highlight_dict', [])
        if isinstance(highlight_dict, str):
            try:
                highlight_dict = json.loads(highlight_dict)
            except (json.JSONDecodeError, TypeError):
                highlight_dict = []

        fragment_core = params.get('fragment_core', '')
        root_title = params.get('root_title', 'SAR Tree')
        max_level = int(params.get('max_level', 5))

        if not fragment_core:
            raise ValueError("fragment_core is required for tree generation")

        # 构造 argparse Namespace 对象
        # tree.py expects Combine_Table_info.csv etc. directly in work_folder,
        # and these files are inside the SAR_Results subdirectory
        work_folder = str(temp_output_dir / "SAR_Results")
        args = Namespace(
            fragment_core=fragment_core,
            rootTitle=root_title,
            workFolder=work_folder,
            inputFile=input_file,
            maxLevel=max_level,
            treeContent=str(tree_content),
            highlightDict=str(highlight_dict),
        )

        await progress_callback.update_progress(20, "Running SAR tree generation")
        logger.info("Task %s: Starting SAR tree generation", task_id)

        # 在线程池中运行
        loop = asyncio.get_event_loop()
        await loop.run_in_executor(self.thread_executor, self._run_tree_sync, args)

        await progress_callback.update_progress(85, "SAR tree generation complete")

    async def submit_task(self, task_id: str, job_dir: str) -> bool:
        """提交新任务"""
        if not self.is_running:
            logger.warning("TaskProcessor is not running, cannot submit task %s", task_id)
            return False

        if task_id in self.active_tasks:
            logger.warning("Task %s is already running", task_id)
            return False

        try:
            task = asyncio.create_task(
                self.process_sarm_task(task_id, job_dir)
            )
            self.active_tasks[task_id] = task

            logger.info("Task %s submitted successfully", task_id)
            return True

        except Exception as e:
            logger.error("Failed to submit task %s: %s", task_id, e)
            return False

    async def cancel_task(self, task_id: str) -> bool:
        """取消任务"""
        if task_id not in self.active_tasks:
            logger.warning("Task %s not found in active tasks", task_id)
            return False

        try:
            task = self.active_tasks[task_id]
            task.cancel()

            try:
                await task
            except asyncio.CancelledError:
                pass

            connection = await self.get_db_connection()
            if connection:
                try:
                    await connection.execute(
                        "UPDATE tasks SET status = $1, finished_at = NOW() WHERE id = $2",
                        "cancelled", task_id
                    )
                finally:
                    await self.release_db_connection(connection)

            del self.active_tasks[task_id]
            logger.info("Task %s cancelled successfully", task_id)
            return True

        except Exception as e:
            logger.error("Failed to cancel task %s: %s", task_id, e)
            return False

    def get_active_tasks(self) -> list:
        """获取活动任务列表"""
        return list(self.active_tasks.keys())

    def get_task_count(self) -> int:
        """获取活动任务数量"""
        return len(self.active_tasks)

    def get_task_progress(self, task_id: str) -> Optional[Dict[str, Any]]:
        """获取特定任务的进度信息"""
        return self.task_progress.get(task_id)

    def get_all_tasks_progress(self) -> Dict[str, Dict[str, Any]]:
        """获取所有任务的进度信息"""
        return self.task_progress.copy()

    async def shutdown(self):
        """关闭任务处理器"""
        logger.info("Shutting down AsyncTaskProcessor...")
        self.is_running = False

        if self.polling_task:
            logger.info("Stopping database polling...")
            self.polling_task.cancel()
            try:
                await self.polling_task
            except asyncio.CancelledError:
                pass

        for task_id, task in self.active_tasks.items():
            logger.info("Cancelling task: %s", task_id)
            task.cancel()

        if self.active_tasks:
            await asyncio.gather(*self.active_tasks.values(), return_exceptions=True)

        if self._db_pool:
            logger.info("Closing database connection pool...")
            await self._db_pool.close()
            self._db_pool = None

        if self.thread_executor:
            logger.info("Shutting down thread executor...")
            self.thread_executor.shutdown(wait=False)

        logger.info("AsyncTaskProcessor shutdown complete")

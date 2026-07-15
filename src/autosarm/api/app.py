"""
FastAPI 应用工厂

创建和配置 FastAPI 应用实例
"""

import logging
import os
from contextlib import asynccontextmanager
from typing import Optional

from fastapi import FastAPI, HTTPException, Request
from fastapi.responses import JSONResponse
from fastapi.exceptions import RequestValidationError
from fastapi.middleware.cors import CORSMiddleware

from autosarm.config.logging import setup_logging
from autosarm.tasks.processor import AsyncTaskProcessor

logger = logging.getLogger(__name__)

# 全局异步任务处理器
_async_processor: Optional[AsyncTaskProcessor] = None


def get_async_processor() -> AsyncTaskProcessor:
    """获取异步任务处理器实例"""
    if _async_processor is None:
        raise RuntimeError("Async processor not initialized")
    return _async_processor


@asynccontextmanager
async def lifespan(app: FastAPI):
    """应用生命周期管理"""
    global _async_processor

    # —— 应用启动时执行 ——
    logger.info("Starting AutoSARM API...")

    # 初始化异步任务处理器
    logger.info("Initializing async task processor...")
    _async_processor = AsyncTaskProcessor()

    # —— 数据库轮询：默认关闭（ADR 0012 P2）——
    #
    # 集群里由 Argo Workflows 领活。**注意 autoSARM 的情况和另外三个 worker 不同**：它在集群里
    # 从来就没有 Deployment、也没有镜像 —— 这个 worker **一次都没跑过**。所以这里没有"回滚"
    # 一说，默认关闭是唯一正确的默认值。
    #
    # 🔴 为什么必须默认关掉：
    # compute-foundry operator 的领活条件是 `status='pending' AND workflow_name IS NULL`。它提交
    # Workflow 后先写 `workflow_name`，`status` 要等下一次 project() 才变成 `processing` ——
    # 中间有个最长 10 秒的窗口，行仍然是 `pending`。任何人在本地把这个应用跑起来并把 DB 指向
    # 生产库，就会在这个窗口里把任务领走 → **同一个任务跑两遍。**
    #
    # 讽刺的是，四个 worker 里 autoSARM 的领活是**唯一写对的**（真事务 + FOR UPDATE SKIP
    # LOCKED）。但它照样要删 —— Argo 领活之后，已经没有东西需要加锁了。
    if os.environ.get("LEGACY_POLLER", "0") == "1":
        logger.warning(
            "LEGACY_POLLER=1 —— 启动数据库轮询。这是 ADR 0012 之前的模式；"
            "而 autoSARM 从未以这种模式在集群里跑过。除非你很清楚自己在做什么，否则不该看到这行。"
        )
        await _async_processor.start_polling()
    else:
        logger.info(
            "数据库轮询已禁用（ADR 0012：调度由 Argo Workflows + compute-foundry operator 接管）。"
            "本进程只提供 HTTP 接口；计算走 `python -m autosarm.steps run`。"
        )

    logger.info("AutoSARM API startup complete")

    yield

    # —— 应用关闭时执行 ——
    logger.info("Shutting down AutoSARM API...")
    if _async_processor:
        await _async_processor.shutdown()
    logger.info("AutoSARM API shutdown complete")


def create_app() -> FastAPI:
    """
    创建 FastAPI 应用实例
    
    Returns:
        配置完成的 FastAPI 应用
    """
    # 设置日志
    setup_logging(level="INFO")

    app = FastAPI(
        lifespan=lifespan,
        title="AutoSARM API",
        description="Automated Structure-Activity Relationship Matrix service. "
                    "Workers poll PostgreSQL for pending tasks and process them automatically.",
        version="1.0.0",
        docs_url="/docs",
        redoc_url="/redoc",
        openapi_url="/openapi.json",
    )

    # 添加 CORS 中间件
    app.add_middleware(
        CORSMiddleware,
        allow_origins=["*"],
        allow_credentials=True,
        allow_methods=["*"],
        allow_headers=["*"],
    )

    # 注册异常处理器
    _register_exception_handlers(app)

    # 注册路由
    _register_routes(app)

    return app


def _register_exception_handlers(app: FastAPI):
    """注册全局异常处理器"""

    @app.exception_handler(RequestValidationError)
    async def validation_exception_handler(request: Request, exc: RequestValidationError):
        logger.error(f"Validation error: {exc}")
        return JSONResponse(
            status_code=422,
            content={"detail": exc.errors(), "body": exc.body}
        )

    @app.exception_handler(HTTPException)
    async def http_exception_handler(request: Request, exc: HTTPException):
        logger.error(f"HTTP error: {exc.detail}")
        return JSONResponse(
            status_code=exc.status_code,
            content={"detail": exc.detail}
        )

    @app.exception_handler(Exception)
    async def general_exception_handler(request: Request, exc: Exception):
        logger.error(f"Unexpected error: {str(exc)}")
        return JSONResponse(
            status_code=500,
            content={"detail": "Internal server error"}
        )


def _register_routes(app: FastAPI):
    """注册路由"""
    from autosarm.api.routes import health

    # 健康检查路由
    app.include_router(health.router, tags=["Health"])

    # 根路由
    @app.get("/")
    async def root():
        """根路径"""
        return {
            "message": "AutoSARM API",
            "version": "1.0.0",
            "docs": "/docs",
        }


# 创建默认应用实例（用于 uvicorn 直接引用）
app = create_app()

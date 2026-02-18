"""
配置加载器 - 从 settings.yaml 加载配置
支持环境变量覆盖: SARM_<SECTION>_<KEY>
"""

import os
from pathlib import Path
from typing import Any, Dict, Optional
from dataclasses import dataclass, field

import yaml

# 尝试加载 .env 文件
try:
    from dotenv import load_dotenv
    load_dotenv()
except ImportError:
    pass


# ============ CPU 核心数自动检测 ============

def _detect_cpu_cores() -> int:
    """
    自动检测可用 CPU 核心数，返回 80% 的核心数（向下取整）
    
    在 Docker 容器中，Python 3.10 的 os.cpu_count() 可能返回宿主机的全部核心数
    而非 Docker 限制的 CPU 数量。因此优先从 cgroup 读取实际限制。
    
    检测优先级：
    1. cgroup v2: /sys/fs/cgroup/cpu.max
    2. cgroup v1: /sys/fs/cgroup/cpu/cpu.cfs_quota_us + cpu.cfs_period_us
    3. os.sched_getaffinity (受 taskset/cpuset 影响)
    4. os.cpu_count() (fallback，可能不受 Docker --cpus 限制)
    """
    try:
        # 尝试 cgroup v2 (现代 Linux / Docker)
        cgroup_v2 = Path("/sys/fs/cgroup/cpu.max")
        if cgroup_v2.exists():
            content = cgroup_v2.read_text().strip()
            parts = content.split()
            if parts[0] != "max":
                quota = int(parts[0])
                period = int(parts[1])
                cgroup_cpus = quota / period
                cores = max(1, int(cgroup_cpus * 0.8))
                return cores

        # 尝试 cgroup v1 (旧版 Docker)
        cgroup_v1_quota = Path("/sys/fs/cgroup/cpu/cpu.cfs_quota_us")
        cgroup_v1_period = Path("/sys/fs/cgroup/cpu/cpu.cfs_period_us")
        if cgroup_v1_quota.exists() and cgroup_v1_period.exists():
            quota = int(cgroup_v1_quota.read_text().strip())
            period = int(cgroup_v1_period.read_text().strip())
            if quota > 0 and period > 0:
                cgroup_cpus = quota / period
                cores = max(1, int(cgroup_cpus * 0.8))
                return cores

        # 尝试 sched_getaffinity (受 cpuset 限制影响)
        try:
            cpu_count = len(os.sched_getaffinity(0))
        except (AttributeError, OSError):
            cpu_count = os.cpu_count() or 4

        # 使用 80% 的 CPU 核心，向下取整，最少 1 核
        cores = max(1, int(cpu_count * 0.8))
        return cores
    except Exception:
        return 4  # 出错时使用安全的默认值


# 应用启动时初始化 CPU 核心数（只计算一次）
DEFAULT_CPU_CORES: int = _detect_cpu_cores()


def get_default_cores() -> int:
    """
    获取默认 CPU 核心数
    
    优先级：
    1. CPU_CORES 环境变量（允许手动覆盖）
    2. 自动检测的 80% CPU 核心数
    """
    env_cores = os.environ.get('CPU_CORES')
    if env_cores:
        try:
            return int(env_cores)
        except ValueError:
            pass
    return DEFAULT_CPU_CORES


# 配置文件搜索路径
def _find_settings_file() -> Optional[Path]:
    """查找配置文件"""
    search_paths = [
        Path(__file__).resolve().parent / "settings.yaml",  # 包内
        Path(__file__).resolve().parent.parent.parent.parent / "config" / "settings.yaml",
        Path.cwd() / "config" / "settings.yaml",
        Path.cwd() / "settings.yaml",
    ]

    for path in search_paths:
        if path.exists():
            return path
    return None


# 缓存配置
_settings_cache: Optional[Dict[str, Any]] = None


def _load_yaml() -> Dict[str, Any]:
    """加载 YAML 配置文件"""
    settings_file = _find_settings_file()
    if settings_file and settings_file.exists():
        with open(settings_file, 'r', encoding='utf-8') as f:
            return yaml.safe_load(f) or {}
    return {}


def _get_env_override(section: str, key: str) -> Optional[str]:
    """获取环境变量覆盖值，前缀: SARM_"""
    env_key = f"SARM_{section.upper()}_{key.upper()}"
    return os.environ.get(env_key)


def _apply_env_overrides(config: Dict[str, Any]) -> Dict[str, Any]:
    """应用环境变量覆盖"""
    for section, values in config.items():
        if isinstance(values, dict):
            for key, value in values.items():
                env_value = _get_env_override(section, key)
                if env_value is not None:
                    if isinstance(value, bool):
                        config[section][key] = env_value.lower() in ('true', '1', 'yes')
                    elif isinstance(value, int):
                        config[section][key] = int(env_value)
                    elif isinstance(value, float):
                        config[section][key] = float(env_value)
                    else:
                        config[section][key] = env_value
    return config


def get_settings() -> Dict[str, Any]:
    """获取完整配置（带缓存）"""
    global _settings_cache
    if _settings_cache is None:
        _settings_cache = _apply_env_overrides(_load_yaml())
    return _settings_cache


def reload_settings() -> Dict[str, Any]:
    """重新加载配置"""
    global _settings_cache
    _settings_cache = None
    return get_settings()


def get(section: str, key: str = None, default: Any = None) -> Any:
    """
    获取配置值
    
    Args:
        section: 配置节名称 (如 'database', 'storage')
        key: 配置键名 (可选，不提供则返回整个节)
        default: 默认值
    """
    cfg = get_settings()
    section_data = cfg.get(section, {})

    if key is None:
        return section_data if section_data else default

    # 支持嵌套键 (如 'pool.min_size')
    if '.' in key:
        keys = key.split('.')
        result = section_data
        for k in keys:
            if isinstance(result, dict):
                result = result.get(k)
            else:
                return default
        return result if result is not None else default

    return section_data.get(key, default)


# ============ 类型安全的配置类 ============

@dataclass
class ServerSettings:
    """服务器配置"""
    host: str = "0.0.0.0"
    port: int = 8030
    title: str = "AutoSARM API"
    description: str = "Automated Structure-Activity Relationship Matrix service"
    version: str = "1.0.0"

    @classmethod
    def from_config(cls) -> "ServerSettings":
        return cls(
            host=get('server', 'host', cls.host),
            port=get('server', 'port', cls.port),
            title=get('server', 'title', cls.title),
            description=get('server', 'description', cls.description),
            version=get('server', 'version', cls.version),
        )


@dataclass
class DatabaseSettings:
    """数据库配置"""
    host: str = "127.0.0.1"
    port: int = 5432
    user: str = "admin"
    password: str = "secret"
    database: str = "mydatabase"
    pool_min_size: int = 1
    pool_max_size: int = 5

    @classmethod
    def from_config(cls) -> "DatabaseSettings":
        return cls(
            host=os.getenv('DB_HOST', get('database', 'host', cls.host)),
            port=int(os.getenv('DB_PORT', get('database', 'port', cls.port))),
            user=os.getenv('DB_USER', get('database', 'user', cls.user)),
            password=os.getenv('DB_PASSWORD', get('database', 'password', cls.password)),
            database=os.getenv('DB_NAME', get('database', 'database', cls.database)),
            pool_min_size=get('database', 'pool.min_size', cls.pool_min_size),
            pool_max_size=get('database', 'pool.max_size', cls.pool_max_size),
        )

    def to_dict(self) -> Dict[str, Any]:
        """转换为字典格式"""
        return {
            "host": self.host,
            "port": self.port,
            "user": self.user,
            "password": self.password,
            "database": self.database,
        }


@dataclass
class StorageSettings:
    """存储配置 (SeaweedFS)"""
    api_type: str = "filer"
    filer_endpoint: str = "http://localhost:8888"
    bucket: str = "astramolecula"
    s3_endpoint: str = "http://localhost:8333"
    access_key: str = ""
    secret_key: str = ""
    temp_dir: str = "/tmp/autosarm"
    presigned_url_expires: int = 3600

    @classmethod
    def from_config(cls) -> "StorageSettings":
        return cls(
            api_type=os.getenv('SEAWEED_API_TYPE', get('storage', 'api_type', cls.api_type)),
            filer_endpoint=os.getenv('SEAWEED_FILER_ENDPOINT', get('storage', 'filer_endpoint', cls.filer_endpoint)),
            bucket=os.getenv('SEAWEED_BUCKET', get('storage', 'bucket', cls.bucket)),
            s3_endpoint=os.getenv('SEAWEED_S3_ENDPOINT', get('storage', 's3_endpoint', cls.s3_endpoint)),
            access_key=os.getenv('SEAWEED_ACCESS_KEY', get('storage', 'access_key', cls.access_key)),
            secret_key=os.getenv('SEAWEED_SECRET_KEY', get('storage', 'secret_key', cls.secret_key)),
            temp_dir=os.getenv('TEMP_DIR', get('storage', 'temp_dir', cls.temp_dir)),
            presigned_url_expires=int(os.getenv('PRESIGNED_URL_EXPIRES', get('storage', 'presigned_url_expires', cls.presigned_url_expires))),
        )

    def get_temp_path(self) -> Path:
        """获取临时目录路径"""
        path = Path(self.temp_dir)
        path.mkdir(parents=True, exist_ok=True)
        return path


@dataclass
class TaskProcessorSettings:
    """任务处理器配置"""
    poll_interval: int = 30

    @classmethod
    def from_config(cls) -> "TaskProcessorSettings":
        return cls(
            poll_interval=int(os.getenv('POLL_INTERVAL', get('task_processor', 'poll_interval', cls.poll_interval))),
        )


@dataclass
class Settings:
    """全局配置"""
    server: ServerSettings = field(default_factory=ServerSettings.from_config)
    database: DatabaseSettings = field(default_factory=DatabaseSettings.from_config)
    storage: StorageSettings = field(default_factory=StorageSettings.from_config)
    task_processor: TaskProcessorSettings = field(default_factory=TaskProcessorSettings.from_config)

    @classmethod
    def load(cls) -> "Settings":
        """加载所有配置"""
        return cls(
            server=ServerSettings.from_config(),
            database=DatabaseSettings.from_config(),
            storage=StorageSettings.from_config(),
            task_processor=TaskProcessorSettings.from_config(),
        )


# 单例配置实例 (延迟加载)
_settings: Optional[Settings] = None


def settings() -> Settings:
    """获取配置单例"""
    global _settings
    if _settings is None:
        _settings = Settings.load()
    return _settings

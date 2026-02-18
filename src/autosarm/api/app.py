"""
FastAPI 应用工厂

创建和配置 FastAPI 应用实例
"""

import logging
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

    # 启动数据库轮询
    logger.info("Starting database polling for sarm_analysis tasks...")
    await _async_processor.start_polling()

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

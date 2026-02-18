"""
健康检查路由
"""

from fastapi import APIRouter, Response

router = APIRouter()


@router.get("/health")
async def health_check():
    """
    健康检查端点
    
    Returns:
        服务状态信息
    """
    return {
        "status": "healthy",
        "service": "autosarm",
        "version": "1.0.0",
    }


@router.get("/health/ready")
async def readiness_check():
    """
    就绪检查端点
    
    检查服务是否准备好接收请求
    """
    return {"status": "ready"}


@router.get("/health/live")
async def liveness_check():
    """
    存活检查端点
    
    用于 Docker healthcheck
    """
    return Response(status_code=200)

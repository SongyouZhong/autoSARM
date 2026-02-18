"""
PostgreSQL 数据库连接模块

支持异步连接池 (asyncpg)
"""

import logging
from typing import Optional

import asyncpg

from autosarm.config.settings import settings

logger = logging.getLogger(__name__)

# 异步连接池
_async_pool: Optional[asyncpg.Pool] = None


async def get_async_pool() -> asyncpg.Pool:
    """
    获取或创建异步连接池
    
    Returns:
        asyncpg 连接池实例
    """
    global _async_pool

    if _async_pool is None:
        try:
            db_settings = settings().database

            logger.info("Creating PostgreSQL async connection pool...")
            _async_pool = await asyncpg.create_pool(
                host=db_settings.host,
                port=db_settings.port,
                user=db_settings.user,
                password=db_settings.password,
                database=db_settings.database,
                min_size=db_settings.pool_min_size,
                max_size=db_settings.pool_max_size,
            )
            logger.info("PostgreSQL async connection pool created successfully")
        except Exception as e:
            logger.error(f"Failed to create PostgreSQL async connection pool: {e}")
            raise

    return _async_pool


async def get_async_connection() -> asyncpg.Connection:
    """
    从异步连接池获取一个连接
    
    Returns:
        asyncpg 连接实例
    """
    pool = await get_async_pool()
    return await pool.acquire()


async def release_async_connection(conn: asyncpg.Connection):
    """
    释放异步连接回连接池
    
    Args:
        conn: 要释放的连接
    """
    pool = await get_async_pool()
    await pool.release(conn)


async def close_pool():
    """关闭连接池"""
    global _async_pool

    if _async_pool is not None:
        logger.info("Closing PostgreSQL connection pool...")
        await _async_pool.close()
        _async_pool = None
        logger.info("PostgreSQL connection pool closed")

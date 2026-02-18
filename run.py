#!/usr/bin/env python3
"""
开发环境快速启动脚本

使用方法:
    python run.py              # 启动 API 服务 (默认 --reload)
    python run.py serve        # 启动 API 服务
    python run.py serve --port 8030  # 指定端口
    python run.py sarm         # 直接运行 SARM 矩阵生成 (CLI 模式)
    python run.py tree         # 直接运行 SAR 树生成 (CLI 模式)
"""

import os
import sys

# 将 src 目录添加到 Python 路径
src_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "src")
if src_path not in sys.path:
    sys.path.insert(0, src_path)


def main():
    """主入口函数"""
    import argparse
    import uvicorn

    parser = argparse.ArgumentParser(
        prog="autosarm",
        description="AutoSARM Service - 开发环境快速启动"
    )

    subparsers = parser.add_subparsers(dest="command", help="可用命令")

    # serve 子命令 - 启动 API 服务（Worker 模式）
    serve_parser = subparsers.add_parser("serve", help="启动 API 服务 (Worker 模式，自动轮询数据库任务)")
    serve_parser.add_argument(
        "--host",
        default="0.0.0.0",
        help="绑定主机 (默认: 0.0.0.0)"
    )
    serve_parser.add_argument(
        "--port",
        type=int,
        default=8030,
        help="绑定端口 (默认: 8030)"
    )
    serve_parser.add_argument(
        "--reload",
        action="store_true",
        default=True,
        help="启用自动重载 (开发模式默认开启)"
    )
    serve_parser.add_argument(
        "--no-reload",
        action="store_true",
        help="禁用自动重载"
    )
    serve_parser.add_argument(
        "--log-level",
        default="info",
        choices=["debug", "info", "warning", "error"],
        help="日志级别 (默认: info)"
    )

    # sarm 子命令 - 直接运行 SARM 矩阵生成
    sarm_parser = subparsers.add_parser("sarm", help="直接运行 SARM 矩阵生成 (CLI 模式)")
    sarm_parser.add_argument("--csvFile", required=True, help="输入 CSV 文件路径")
    sarm_parser.add_argument("--type", default="smiles", choices=["smiles", "scaffold"])
    sarm_parser.add_argument("--column", nargs="+", required=True, help="活性列名称")
    sarm_parser.add_argument("--log", type=int, default=0)
    sarm_parser.add_argument("--minimumSite1", type=float, default=3)
    sarm_parser.add_argument("--minimumSite2", type=float, default=3)
    sarm_parser.add_argument("--n_jobs", type=int, default=8)
    sarm_parser.add_argument("--save_folder", default="SAR_Results")
    sarm_parser.add_argument("--csv2excel", type=int, default=0)

    # tree 子命令 - 直接运行 SAR 树生成
    tree_parser = subparsers.add_parser("tree", help="直接运行 SAR 树生成 (CLI 模式)")
    tree_parser.add_argument("--fragment_core", required=True)
    tree_parser.add_argument("--rootTitle", required=True)
    tree_parser.add_argument("--workFolder", required=True)
    tree_parser.add_argument("--inputFile", default="input.csv")
    tree_parser.add_argument("--maxLevel", type=int, default=5)
    tree_parser.add_argument("--treeContent", type=str, default="['double-cut']")
    tree_parser.add_argument("--highlightDict", type=str, default="[]")

    args = parser.parse_args()

    if args.command is None:
        # 默认启动 API 服务
        print("Starting AutoSARM API server (dev mode, reload enabled)")
        print(f"  API docs: http://localhost:8030/docs")
        print(f"  Health:   http://localhost:8030/health")
        print(f"  Press Ctrl+C to stop\n")
        uvicorn.run(
            "autosarm.api.app:create_app",
            factory=True,
            host="0.0.0.0",
            port=8030,
            reload=True,
            log_level="info",
        )
    elif args.command == "serve":
        reload_enabled = args.reload and not args.no_reload
        print(f"Starting AutoSARM API server")
        print(f"  Host: {args.host}")
        print(f"  Port: {args.port}")
        print(f"  Reload: {'on' if reload_enabled else 'off'}")
        print(f"  API docs: http://{args.host if args.host != '0.0.0.0' else 'localhost'}:{args.port}/docs")
        print(f"  Press Ctrl+C to stop\n")
        uvicorn.run(
            "autosarm.api.app:create_app",
            factory=True,
            host=args.host,
            port=args.port,
            reload=reload_enabled,
            workers=1,
            log_level=args.log_level,
        )
    elif args.command == "sarm":
        from autosarm.cli.create_sarm import run_sarm
        run_sarm(args)
    elif args.command == "tree":
        from autosarm.cli.create_tree import run_tree
        run_tree(args)


if __name__ == "__main__":
    main()

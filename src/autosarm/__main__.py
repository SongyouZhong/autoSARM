"""
AutoSARM 包入口

支持:
    python -m autosarm serve           # 启动 Worker 服务
    python -m autosarm sarm --help     # SARM CLI
    python -m autosarm tree --help     # Tree CLI
"""

import sys
import os

# 确保 src 目录在路径中
src_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..")
if src_path not in sys.path:
    sys.path.insert(0, src_path)


def main():
    """包级入口"""
    args = sys.argv[1:]

    if not args or args[0] == "serve":
        # Worker 服务模式
        import argparse
        import uvicorn

        parser = argparse.ArgumentParser(prog="autosarm serve")
        parser.add_argument("--host", default="0.0.0.0")
        parser.add_argument("--port", type=int, default=8030)
        parser.add_argument("--reload", action="store_true")
        parser.add_argument("--log-level", default="info")

        # 移除 'serve' 如果存在
        serve_args = [a for a in args if a != "serve"]
        parsed = parser.parse_args(serve_args)

        uvicorn.run(
            "autosarm.api.app:create_app",
            factory=True,
            host=parsed.host,
            port=parsed.port,
            reload=parsed.reload,
            workers=1,
            log_level=parsed.log_level,
        )
    else:
        # 委托给原始 CLI
        from autosarm.cli.main import main as cli_main
        cli_main()


if __name__ == "__main__":
    main()

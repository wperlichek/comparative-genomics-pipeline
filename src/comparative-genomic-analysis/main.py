import asyncio
from .config import logging_config


async def main():
    logging_config.setup_logging()
    pass


if __name__ == "__main__":
    asyncio.run(main())

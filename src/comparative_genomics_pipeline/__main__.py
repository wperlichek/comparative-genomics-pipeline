import asyncio, logging
from .config import logging_config

logger = logging.getLogger(__name__)

async def async_main():
    pass


def main():
    logging_config.setup_logging()
    logger.info("Starting pipeline...")
    asyncio.run(async_main())


if __name__ == "__main__":
    main()

from logging import DEBUG, INFO, WARNING, ERROR, CRITICAL, basicConfig, getLogger

basicConfig(encoding='utf-8')
logger = getLogger(__name__)

## default setting (change it)
logger.setLevel(INFO)
#logger.setLevel(DEBUG)

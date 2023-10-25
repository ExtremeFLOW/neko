import logging
# Retrieved from https://stackoverflow.com/questions/384076/how-can-i-color-python-logging-output

class CustomFormatter(logging.Formatter):

    grey = "\x1b[38;20m"
    yellow = "\x1b[33;20m"
    red = "\x1b[31;20m"
    bold_red = "\x1b[31;1m"
    reset = "\x1b[0m"
    format = "%(asctime)s - %(name)s - %(levelname)s - %(message)s (%(filename)s:%(lineno)d)"

    FORMATS = {
        logging.DEBUG: grey + format + reset,
        logging.INFO: grey + format + reset,
        logging.WARNING: yellow + format + reset,
        logging.ERROR: red + format + reset,
        logging.CRITICAL: bold_red + format + reset
    }

    def format(self, record):
        log_fmt = self.FORMATS.get(record.levelno)
        formatter = logging.Formatter(log_fmt)
        return formatter.format(record)


class logger_c():
    def __init__(self, level: None, comm: None):

        self.level = level
        self.comm = comm

        #Instanciate
        logger = logging.getLogger(__name__)
        logger.setLevel(level)

        # create console handler with a higher log level
        ch = logging.StreamHandler()
        ch.setLevel(level)

        ch.setFormatter(CustomFormatter())

        logger.addHandler(ch)

        self.log = logger

    def write(self,level, message):
        
        comm = self.comm
        rank     = comm.Get_rank()
        size     = comm.Get_size()

        if level == "debug":
            if rank == 0: self.log.debug(message)
        
        if level == "info":
            if rank == 0: self.log.info(message)

        if level == "warning":
            if rank == 0: self.log.warning(message)

        if level == "error":
            if rank == 0: self.log.error(message)

        if level == "critical":
            if rank == 0: self.log.critical(message)


class CookieCutter:
    """
    Cookie Cutter wrapper
    """

    @staticmethod
    def get_default_mem_mb() -> int:
        return int("1024")

    @staticmethod
    def get_log_dir() -> str:
        return "logs/cluster"

    @staticmethod
    def get_default_queue() -> str:
        return "medium"

    @staticmethod
    def get_default_project() -> str:
        return "re_gecip_cancer_pan"

    @staticmethod
    def get_lsf_unit_for_limits() -> str:
        return "MB"

    @staticmethod
    def get_unknwn_behaviour() -> str:
        return "kill"

    @staticmethod
    def get_zombi_behaviour() -> str:
        return "kill"

    @staticmethod
    def get_latency_wait() -> float:
        return float("60")

    @staticmethod
    def get_wait_between_tries() -> float:
        return float("0.001")

    @staticmethod
    def get_max_status_checks() -> int:
        return int("1")

    @staticmethod
    def jobscript_timeout() -> int:
        return int("10")

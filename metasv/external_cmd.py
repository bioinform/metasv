import time
import shlex
import subprocess
from threading import Timer
import unittest
import logging

class TimedExternalCmd:
    def __init__(self, cmd, logger):
        self.cmd = shlex.split(cmd)
        self.p = None
        self.did_timeout = False
        self.logger = logger
    def enforce_timeout(self):
        self.p.terminate()
        self.did_timeout = True
    def run(self, cmd_log_fd_out=None, cmd_log_fd_err=None, timeout=None):
        self.logger.info("Running %s with arguments %s" % (self.cmd[0].upper(), str(self.cmd[1::])))
        cmd_log_fd_err = cmd_log_fd_err or cmd_log_fd_out
        self.p = subprocess.Popen(self.cmd, stderr=cmd_log_fd_err, stdout=cmd_log_fd_out)
        start_time = time.time()
        if timeout:
            t = Timer(timeout, self.enforce_timeout)
            t.start()
        self.p.wait()
        if timeout:
            t.cancel()
            if self.did_timeout:
                self.logger.error("Timed out after %d seconds", timeout)
                return None
        retcode = self.p.returncode
        self.logger.info("Returned code %d (%g seconds)" % (retcode, time.time() - start_time))
        return retcode


class TestTimedExternalCmd(unittest.TestCase):
    def test_run_complete(self):
        cmd = TimedExternalCmd("sleep 1", self.logger)
        self.assertEqual(cmd.run(timeout = 2), 0)
        self.assertFalse(cmd.did_timeout)
        return

    def test_run_timeout(self):
        start_tick = time.time()
        cmd = TimedExternalCmd("sleep 2", self.logger)
        cmd.run(timeout = 1)
        run_time = time.time() - start_tick
        self.assertTrue(cmd.did_timeout)
        self.assertAlmostEqual(run_time, 1, delta=0.2)
        return

    def test_run_no_timeout(self):
        cmd = TimedExternalCmd("sleep 1", self.logger)
        retcode = cmd.run()
        self.assertEqual(cmd.run(), 0)
        self.assertFalse(cmd.did_timeout)
        return

    def test_run_fail(self):
        cmd = TimedExternalCmd("sleep 1 2 3", self.logger)
        retcode = cmd.run(timeout = 1)
        self.assertIsNotNone(retcode)
        self.assertIsNot(retcode, 0)
        return

    logger = None


if __name__ == '__main__':
    FORMAT = '%(levelname)s %(asctime)-15s %(name)-20s %(message)s'
    logging.basicConfig(level=logging.INFO, format=FORMAT)
    TestTimedExternalCmd.logger = logging.getLogger(__name__)
    unittest.main()

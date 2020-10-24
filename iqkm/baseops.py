#!/usr/local/bin/python3

import os
import logging


class file:
    @staticmethod
    def isfile(f):
        return os.path.exists(f)

    @staticmethod
    def exists(f):
        return os.path.exists(f)

    @staticmethod
    def isdir(d, create=True):
        """
        check if dir exists and create if not
        """
        if not os.path.isdir(d):
            if create:
                try:
                    os.makedirs(d)
                    return True
                except OSError as e:
                    logging.warning(f"Could not create dir: {d}\n{e}")
                    return False
            else:
                return False
        else:
            return True

    @staticmethod
    def isnewer(fileA, fileB):
        """
        Check if file A is newer than file B
        """
        if not file.isfile(fileA):
            logging.warning("{} is not a file".format(fileA))
            return False

        if not file.isfile(fileB):
            # return true bc fileB does not exists!
            return True
        return os.stat(fileA).st_mtime > os.stat(fileB).st_mtime

    @staticmethod
    def touch(path):
        with open(path, "a"):
            os.utime(path, None)

    @staticmethod
    def which(program):
        """
        test if w programm is avaliable

        :param programm: name of an executable

        :rtype: bool
        """
        # taken from
        # https://stackoverflow.com/questions/377017/\
        # test-if-executable-exists-in-python
        def is_exe(fpath):
            return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

        fpath, fname = os.path.split(program)
        if fpath:
            if is_exe(program):
                return program
        else:
            for path in os.environ["PATH"].split(os.pathsep):
                exe_file = os.path.join(path, program)
                if is_exe(exe_file):
                    return exe_file

        return None


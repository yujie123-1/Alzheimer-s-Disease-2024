# coding:utf-8
# yujie
# MAIL : yj.yang1@siat.ac.cn
# DATE : 11/24/2022 4:18 PM

# encoding: utf-8
# module time
# from (built-in)
# by generator 1.147
"""
This module provides various functions to manipulate time values.

There are two standard representations of time.  One is the number
of seconds since the Epoch, in UTC (a.k.a. GMT).  It may be an integer
or a floating point number (to represent fractions of seconds).
The Epoch is system-defined; on Unix, it is generally January 1st, 1970.
The actual value can be retrieved by calling gmtime(0).

The other representation is a tuple of 9 integers giving local time.
The tuple items are:
  year (including century, e.g. 1998)
  month (1-12)
  day (1-31)
  hours (0-23)
  minutes (0-59)
  seconds (0-59)
  weekday (0-6, Monday is 0)
  Julian day (day in the year, 1-366)
  DST (Daylight Savings Time) flag (-1, 0 or 1)
If the DST flag is 0, the time is given in the regular time zone;
if it is 1, the time is given in the DST time zone;
if it is -1, mktime() should guess based on the date and time.
"""
# no imports

# Variables with simple values

altzone = -14400

daylight = 0

timezone = -10800

_STRUCT_TM_ITEMS = 11


# functions

def asctime(p_tuple=None):  # real signature unknown; restored from __doc__
    """
    asctime([tuple]) -> string

    Convert a time tuple to a string, e.g. 'Sat Jun 06 16:26:11 1998'.
    When the time tuple is not present, current time as returned by localtime()
    is used.
    """
    return ""


def ctime(seconds=None):  # known case of time.ctime
    """
    ctime(seconds) -> string

    Convert a time in seconds since the Epoch to a string in local time.
    This is equivalent to asctime(localtime(seconds)). When the time tuple is
    not present, current time as returned by localtime() is used.
    """
    return ""


def get_clock_info(name):  # real signature unknown; restored from __doc__
    """
    get_clock_info(name: str) -> dict

    Get information of the specified clock.
    """
    return {}


def gmtime(seconds=None):  # real signature unknown; restored from __doc__
    """
    gmtime([seconds]) -> (tm_year, tm_mon, tm_mday, tm_hour, tm_min,
                           tm_sec, tm_wday, tm_yday, tm_isdst)

    Convert seconds since the Epoch to a time tuple expressing UTC (a.k.a.
    GMT).  When 'seconds' is not passed in, convert the current time instead.

    If the platform supports the tm_gmtoff and tm_zone, they are available as
    attributes only.
    """
    pass


def localtime(seconds=None):  # real signature unknown; restored from __doc__
    """
    localtime([seconds]) -> (tm_year,tm_mon,tm_mday,tm_hour,tm_min,
                              tm_sec,tm_wday,tm_yday,tm_isdst)

    Convert seconds since the Epoch to a time tuple expressing local time.
    When 'seconds' is not passed in, convert the current time instead.
    """
    pass


def mktime(p_tuple):  # real signature unknown; restored from __doc__
    """
    mktime(tuple) -> floating point number

    Convert a time tuple in local time to seconds since the Epoch.
    Note that mktime(gmtime(0)) will not generally return zero for most
    time zones; instead the returned value will either be equal to that
    of the timezone or altzone attributes on the time module.
    """
    return 0.0


def monotonic():  # real signature unknown; restored from __doc__
    """
    monotonic() -> float

    Monotonic clock, cannot go backward.
    """
    return 0.0


def monotonic_ns():  # real signature unknown; restored from __doc__
    """
    monotonic_ns() -> int

    Monotonic clock, cannot go backward, as nanoseconds.
    """
    return 0


def perf_counter():  # real signature unknown; restored from __doc__
    """
    perf_counter() -> float

    Performance counter for benchmarking.
    """
    return 0.0


def perf_counter_ns():  # real signature unknown; restored from __doc__
    """
    perf_counter_ns() -> int

    Performance counter for benchmarking as nanoseconds.
    """
    return 0


def process_time():  # real signature unknown; restored from __doc__
    """
    process_time() -> float

    Process time for profiling: sum of the kernel and user-space CPU time.
    """
    return 0.0


def process_time_ns(*args, **kwargs):  # real signature unknown
    """
    process_time() -> int

    Process time for profiling as nanoseconds:
    sum of the kernel and user-space CPU time.
    """
    pass


def sleep(seconds):  # real signature unknown; restored from __doc__
    """
    sleep(seconds)

    Delay execution for a given number of seconds.  The argument may be
    a floating point number for subsecond precision.
    """
    pass


def strftime(format, p_tuple=None):  # real signature unknown; restored from __doc__
    """
    strftime(format[, tuple]) -> string

    Convert a time tuple to a string according to a format specification.
    See the library reference manual for formatting codes. When the time tuple
    is not present, current time as returned by localtime() is used.

    Commonly used format codes:

    %Y  Year with century as a decimal number.
    %m  Month as a decimal number [01,12].
    %d  Day of the month as a decimal number [01,31].
    %H  Hour (24-hour clock) as a decimal number [00,23].
    %M  Minute as a decimal number [00,59].
    %S  Second as a decimal number [00,61].
    %z  Time zone offset from UTC.
    %a  Locale's abbreviated weekday name.
    %A  Locale's full weekday name.
    %b  Locale's abbreviated month name.
    %B  Locale's full month name.
    %c  Locale's appropriate date and time representation.
    %I  Hour (12-hour clock) as a decimal number [01,12].
    %p  Locale's equivalent of either AM or PM.

    Other codes may be available on your platform.  See documentation for
    the C library strftime function.
    """
    return ""


def strptime(string, format):  # real signature unknown; restored from __doc__
    """
    strptime(string, format) -> struct_time

    Parse a string to a time tuple according to a format specification.
    See the library reference manual for formatting codes (same as
    strftime()).

    Commonly used format codes:

    %Y  Year with century as a decimal number.
    %m  Month as a decimal number [01,12].
    %d  Day of the month as a decimal number [01,31].
    %H  Hour (24-hour clock) as a decimal number [00,23].
    %M  Minute as a decimal number [00,59].
    %S  Second as a decimal number [00,61].
    %z  Time zone offset from UTC.
    %a  Locale's abbreviated weekday name.
    %A  Locale's full weekday name.
    %b  Locale's abbreviated month name.
    %B  Locale's full month name.
    %c  Locale's appropriate date and time representation.
    %I  Hour (12-hour clock) as a decimal number [01,12].
    %p  Locale's equivalent of either AM or PM.

    Other codes may be available on your platform.  See documentation for
    the C library strftime function.
    """
    return struct_time


def thread_time():  # real signature unknown; restored from __doc__
    """
    thread_time() -> float

    Thread time for profiling: sum of the kernel and user-space CPU time.
    """
    return 0.0


def thread_time_ns(*args, **kwargs):  # real signature unknown
    """
    thread_time() -> int

    Thread time for profiling as nanoseconds:
    sum of the kernel and user-space CPU time.
    """
    pass


def time():  # real signature unknown; restored from __doc__
    """
    time() -> floating point number

    Return the current time in seconds since the Epoch.
    Fractions of a second may be present if the system clock provides them.
    """
    return 0.0


def time_ns():  # real signature unknown; restored from __doc__
    """
    time_ns() -> int

    Return the current time in nanoseconds since the Epoch.
    """
    return 0


# classes

class struct_time(tuple):
    """
    The time value as returned by gmtime(), localtime(), and strptime(), and
     accepted by asctime(), mktime() and strftime().  May be considered as a
     sequence of 9 integers.

     Note that several fields' values are not the same as those defined by
     the C language standard for struct tm.  For example, the value of the
     field tm_year is the actual year, not year - 1900.  See individual
     fields' descriptions for details.
    """

    def __init__(self, *args, **kwargs):  # real signature unknown
        pass

    @staticmethod  # known case of __new__
    def __new__(*args, **kwargs):  # real signature unknown
        """ Create and return a new object.  See help(type) for accurate signature. """
        pass

    def __reduce__(self, *args, **kwargs):  # real signature unknown
        pass

    def __repr__(self, *args, **kwargs):  # real signature unknown
        """ Return repr(self). """
        pass

    tm_gmtoff = property(lambda self: object(), lambda self, v: None, lambda self: None)  # default
    """offset from UTC in seconds"""

    tm_hour = property(lambda self: object(), lambda self, v: None, lambda self: None)  # default
    """hours, range [0, 23]"""

    tm_isdst = property(lambda self: object(), lambda self, v: None, lambda self: None)  # default
    """1 if summer time is in effect, 0 if not, and -1 if unknown"""

    tm_mday = property(lambda self: object(), lambda self, v: None, lambda self: None)  # default
    """day of month, range [1, 31]"""

    tm_min = property(lambda self: object(), lambda self, v: None, lambda self: None)  # default
    """minutes, range [0, 59]"""

    tm_mon = property(lambda self: object(), lambda self, v: None, lambda self: None)  # default
    """month of year, range [1, 12]"""

    tm_sec = property(lambda self: object(), lambda self, v: None, lambda self: None)  # default
    """seconds, range [0, 61])"""

    tm_wday = property(lambda self: object(), lambda self, v: None, lambda self: None)  # default
    """day of week, range [0, 6], Monday is 0"""

    tm_yday = property(lambda self: object(), lambda self, v: None, lambda self: None)  # default
    """day of year, range [1, 366]"""

    tm_year = property(lambda self: object(), lambda self, v: None, lambda self: None)  # default
    """year, for example, 1993"""

    tm_zone = property(lambda self: object(), lambda self, v: None, lambda self: None)  # default
    """abbreviation of timezone name"""

    n_fields = 11
    n_sequence_fields = 9
    n_unnamed_fields = 0
    __match_args__ = (
        'tm_year',
        'tm_mon',
        'tm_mday',
        'tm_hour',
        'tm_min',
        'tm_sec',
        'tm_wday',
        'tm_yday',
        'tm_isdst',
    )


class __loader__(object):
    """
    Meta path import for built-in modules.

        All methods are either class or static methods to avoid the need to
        instantiate the class.
    """

    def create_module(spec):  # reliably restored by inspect
        """ Create a built-in module """
        pass

    def exec_module(module):  # reliably restored by inspect
        """ Exec a built-in module """
        pass

    @classmethod
    def find_module(cls, *args, **kwargs):  # real signature unknown
        """
        Find the built-in module.

                If 'path' is ever specified then the search is considered a failure.

                This method is deprecated.  Use find_spec() instead.
        """
        pass

    @classmethod
    def find_spec(cls, *args, **kwargs):  # real signature unknown
        pass

    @classmethod
    def get_code(cls, *args, **kwargs):  # real signature unknown
        """ Return None as built-in modules do not have code objects. """
        pass

    @classmethod
    def get_source(cls, *args, **kwargs):  # real signature unknown
        """ Return None as built-in modules do not have source code. """
        pass

    @classmethod
    def is_package(cls, *args, **kwargs):  # real signature unknown
        """ Return False as built-in modules are never packages. """
        pass

    @classmethod
    def load_module(cls, *args, **kwargs):  # real signature unknown
        """
        Load the specified module into sys.modules and return it.

            This method is deprecated.  Use loader.exec_module() instead.
        """
        pass

    def module_repr(module):  # reliably restored by inspect
        """
        Return repr for the module.

                The method is deprecated.  The import machinery does the job itself.
        """
        pass

    def __init__(self, *args, **kwargs):  # real signature unknown
        pass

    __weakref__ = property(lambda self: object(), lambda self, v: None, lambda self: None)  # default
    """list of weak references to the object (if defined)"""

    _ORIGIN = 'built-in'
    __dict__ = None  # (!) real value is "mappingproxy({'__module__': '_frozen_importlib', '__doc__': 'Meta path import for built-in modules.\\n\\n    All methods are either class or static methods to avoid the need to\\n    instantiate the class.\\n\\n    ', '_ORIGIN': 'built-in', 'module_repr': <staticmethod(<function BuiltinImporter.module_repr at 0x000001DA29A92320>)>, 'find_spec': <classmethod(<function BuiltinImporter.find_spec at 0x000001DA29A923B0>)>, 'find_module': <classmethod(<function BuiltinImporter.find_module at 0x000001DA29A92440>)>, 'create_module': <staticmethod(<function BuiltinImporter.create_module at 0x000001DA29A924D0>)>, 'exec_module': <staticmethod(<function BuiltinImporter.exec_module at 0x000001DA29A92560>)>, 'get_code': <classmethod(<function BuiltinImporter.get_code at 0x000001DA29A92680>)>, 'get_source': <classmethod(<function BuiltinImporter.get_source at 0x000001DA29A927A0>)>, 'is_package': <classmethod(<function BuiltinImporter.is_package at 0x000001DA29A928C0>)>, 'load_module': <classmethod(<function _load_module_shim at 0x000001DA29A917E0>)>, '__dict__': <attribute '__dict__' of 'BuiltinImporter' objects>, '__weakref__': <attribute '__weakref__' of 'BuiltinImporter' objects>})"


# variables with complex values

tzname = (
    'Russia TZ 2 Standard Time',
    'Russia TZ 2 Daylight Time',
)

__spec__ = None  # (!) real value is "ModuleSpec(name='time', loader=<class '_frozen_importlib.BuiltinImporter'>, origin='built-in')"



import time
import os
start = time.perf_counter()
print("Execute Program")

# list all files in current folder
a = r"E:\AD_Patient\APE\RawData\rs429358"
cwd = os.listdir(r'E:\AD_Patient\APE\RawData\rs429358')
for x in cwd:
    files = a + "\\" + x   # construct path
# change path
os.chdir('E:/AD_Patient/APE/RawData/rs429358')
for file in cwd:
    text = open(file, "r", encoding='utf-8')
    for line in text:
        linesplit = line.split("\t")
        chrom = linesplit[0]
        position = linesplit[1]
        ref = linesplit[2]
        A1 = str(linesplit[4].count("A") + linesplit[4].count("a"))
        T1 = str(linesplit[4].count("T") + linesplit[4].count("t"))
        C1 = str(linesplit[4].count("C") + linesplit[4].count("c"))
        G1 = str(linesplit[4].count("G") + linesplit[4].count("g"))
        APOE1 = "rs429358"
        T1edit = "0"
        C1edit = "0"
        AP0E1edit = "rs429358"
        # count the number of bases per position
        if ref == "A":
            A1 = str(linesplit[4].count(".") + linesplit[4].count(","))
        elif ref == "T":
            T1 = str(linesplit[4].count(".") + linesplit[4].count(","))
        elif ref == "C":
            C1 = str(linesplit[4].count(".") + linesplit[4].count(","))
        elif ref == "G":
            G1 = str(linesplit[4].count(".") + linesplit[4].count(","))
        else:
            continue
        # identify the number of APOE4 of each sample in the raw data
        if int(T1) != 0 and int(C1) == 0:
            APOE1 = "noE4"
        elif int(T1) == 0 and int(C1) != 0:
            APOE1 = "E4/4"
        elif int(T1) != 0 and int(C1) != 0:
            APOE1 = "E4carrier"
        else:
            APOE1 = "NoDetect"
        # correct the data
        if int(C1) == 0:
            T1edit = T1
            C1edit = C1
        elif int(T1) / int(C1) >= 5 and int(C1) != 0:
            T1edit = T1
            C1edit = "0"
        elif int(T1) / int(C1) <= 0.2 and int(C1) != 0:
            T1edit = "0"
            C1edit = C1
        else:
            T1edit = T1
            C1edit = C1
        # identify the number of APOE4 of each sample after correction data
        if int(T1edit) != 0 and int(C1edit) == 0:
            AP0E1edit = "noE4"
        elif int(T1edit) == 0 and int(C1edit) != 0:
            AP0E1edit = "E4/4"
        elif int(T1edit) != 0 and int(C1edit) != 0:
            AP0E1edit = "E4carrier"
        else:
            AP0E1edit = "NoDetect"
        # construct the results into a list
        output1 = "\t".join([file[0:10], chrom, position, A1, T1, C1, G1, APOE1, C1edit, T1edit, AP0E1edit])
        # write the results into a file
        # eg: E:/AD_Patient/APE/rs429358.txt
        fp = open("E:/AD_Patient/APE/rs429358.txt", 'a+')
        print(output1, file=fp)
        fp.close()
# change the files path and check the number of lines in file
# eg:os.chdir('E:/AD_Patient/APE')
os.chdir('E:/AD_Patient/APE')
# check the number of lines in file
file1 = open("rs429358.txt", "r")
file1_lists = file1.readlines()
# print("--------Display the First 5 lines of rs429358-----------")
print("Sample Number Collected at rs429358 Site:" + str(len(file1_lists)))
# check the first few lines in a file
file1 = open("rs429358.txt", "r")
# for i in range(int(5)):
#     a = file1.readline().strip()
#     a = a.split("\t")
#     Rs429358 = a[0:7]
#     print(Rs429358)
# print("", end='\n')
# calculate the second site
# identify the number of APOE4 of each sample
# the second site is not used, the rs429358 SNP site is mainly used,
# but in order to identify APOE genotype, the rs7412 SNP site still is extracted
a = r"E:\AD_Patient\APE\RawData\rs7412"
cwd = os.listdir(r'E:\AD_Patient\APE\RawData\rs7412')
for x in cwd:
    files = a + "\\" + x  # construct path
# change path
os.chdir('E:/AD_Patient/APE/RawData/rs7412')
for file in cwd:
    text = open(file, "r", encoding='utf-8')
    for line in text:
        linesplit = line.split("\t")
        chrom = linesplit[0]
        position = linesplit[1]
        ref = linesplit[2]
        A2 = str(linesplit[4].count("A") + linesplit[4].count("a"))
        T2 = str(linesplit[4].count("T") + linesplit[4].count("t"))
        C2 = str(linesplit[4].count("C") + linesplit[4].count("c"))
        G2 = str(linesplit[4].count("G") + linesplit[4].count("g"))
        T1edit = "0"
        C1edit = "0"
        # count the number of bases per position
        if ref == "A":
            A2 = str(linesplit[4].count(".") + linesplit[4].count(","))
        elif ref == "T":
            T2 = str(linesplit[4].count(".") + linesplit[4].count(","))
        elif ref == "C":
            C2 = str(linesplit[4].count(".") + linesplit[4].count(","))
        elif ref == "G":
            G2 = str(linesplit[4].count(".") + linesplit[4].count(","))
        else:
            continue
        # correct the data
        if int(C2) == 0:
            T2edit = T2
            C2edit = C2
        elif int(T2) / int(C2) >= 5 and int(C2) != 0:
            T2edit = T2
            C2edit = "0"
        elif int(T2) / int(C2) <= 0.2 and int(C2) != 0:
            T2edit = "0"
            C2edit = C2
        else:
            T2edit = T2
            C2edit = C2
        # construct the result into a list
        output2 = "\t".join([file[0:10], chrom, position, A2, T2, C2, G2, C2edit, T2edit])
        # write the result into a file
        fp = open('E:/AD_Patient/APE/rs7412.txt', 'a+')
        print(output2, file=fp)
        fp.close()
# change path
os.chdir('E:/AD_Patient/APE')
# check the number of lines in a file
file2 = open("rs7412.txt", "r")
file2_lists = file2.readlines()
# print("--------Display the Content of rs7412 SNP Site-----------")
print("Sample Number Collected at rs7412 Site:" + str(len(file1_lists)))
# check the first few lines in file
file2 = open("rs7412.txt", "r")
# for i in range(int(5)):
#     a = file1.readline().strip()
#     a = a.split("\t")
#     Rs7412 = a[0:7]
#     print(Rs7412)
print("", end='\n')
# rs429358 SNP file and rs7412 SNP file are reloaded
# identify APOE genotype of each sample
# change path
os.chdir('E:/AD_Patient/APE')
# define the number of lines you want to check
n = int(5)
Rs429358 = open("rs429358.txt", "r", encoding='utf-8')  # open the rs429358 SNP site file
# check the content of the rs429358 SNP site file
print("-----------Show Predict the Number of APOE4 in Each Sample-----------")
print("Display Prediction Result:")
for i in range(n):
    print(Rs429358.readline().strip())
print("", end='\n')
# check sample name of rs429358 and rs7412 site file
# make preparation for subsequent data merge
print("-----------Check Sample Names-----------")
file1 = open("rs429358.txt", "r")
file2 = open("rs7412.txt", "r")
file1_lists = file1.readlines()
file2_lists = file2.readlines()
file3_list = []
file4_list = []
for i in file1_lists:
    temp_list = i.split()
    file3_list.append(str(temp_list[0]))
# print the sample name of rs429358 SNP site in a list
# only the result of the last loop is printed when print is written outside the loop
print("rs429358 Site Sample Names:")
print(file3_list)
for i in file2_lists:
    temp_list = i.split()
    file4_list.append(str(temp_list[0]))
# print the sample name of rs7412 SNP site in a list
# only the result of the last loop is printed when print is written outside the loop
print("rs7412 Site Sample Names:")
print(file4_list)
# print the comparison result between the two SNP sites files
if file3_list == file4_list:
    print("rs429358 and rs7412 site samples were in the same order")
    print("", end='\n')
else:
    print("rs429358 and rs7412 site samples were not in the same order")
    print("", end='\n')
# compare the inconsistent samples in the two site files
print("-----------Check Sample Names in Both Sites-----------")
# print the samples that are in rs429358 but not in rs7412
diff1 = set(file3_list).difference(set(file4_list))
print("samples that are in rs429358 file but not in rs7412 file:" + str(diff1))
# print the samples that are in rs7412 but not in rs429358
diff2 = set(file4_list).difference(set(file3_list))
print("samples that are in rs7412 file but not in rs429358 file:" + str(diff2))
# print the common samples that are in both rs7412 and rs429358
common = set(file3_list).intersection(set(file4_list))
# print(common)
# count sample number
print("Common Samples Number in Both Files:" + str(len(common)))
print("", end='\n')
print("-----------Extract Common Samples in Both Files------------")
file1 = open("rs429358.txt", "r", encoding='utf-8')
for line in file1:
    linesplit = line.split("\t")
    Sample = linesplit[0]
    if Sample in common:
        output1 = "\t".join(linesplit)
        # remove the newline character of output1
        output1 = output1.replace('\n', '').replace('\r', '')
        # write the result into a file
        fp = open('E:/AD_Patient/APE/rs429358Common.txt', 'a+')
        print(output1, file=fp)
        fp.close()
file2 = open("rs7412.txt", "r", encoding='utf-8')
for line in file2:
    linesplit = line.split("\t")
    Sample = linesplit[0]
    if Sample in common:
        output2 = "\t".join(linesplit)
        # remove the newline character of output1
        output2 = output2.replace('\n', '').replace('\r', '')
        # write the result into a file
        fp = open('E:/AD_Patient/APE/rs7412Common.txt', 'a+')
        print(output2, file=fp)
        fp.close()
# check sample names
print("-----------Recheck Sample Names-----------")
file1 = open("rs429358Common.txt", "r")
file2 = open("rs7412Common.txt", "r")
file1_lists = file1.readlines()
file2_lists = file2.readlines()
file3_list = []
file4_list = []
# rs429358 site file
for i in file1_lists:
    temp_list = i.split()
    file3_list.append(str(temp_list[0]))
# print sample of rs429358 site file
# only the result of the last loop is printed when print is written outside the loop
print("Number of Common Samples in rs429358 Site:" + str(len(file3_list)))
print("Common Samples in rs429358 Site:")
print(file3_list)
# rs7412 site file
for i in file2_lists:
    temp_list = i.split()
    file4_list.append(str(temp_list[0]))
# print sample of rs7412 site file
# only the result of the last loop is printed when print is written outside the loop
print("Number of Common Samples in rs7412 Site:" + str(len(file4_list)))
print("Common Samples in rs7412 Site:")
print(file4_list)
# check sample order
# print the result of sample comparison
if file3_list == file4_list:
    print("rs429358 and rs7412 site samples were in the same order")
    print("", end='\n')
else:
    print("rs429358 and rs7412 site samples were not in the same order")
    print("", end='\n')
# merge two files in the same order
# read the rs429358 site file to concatenate
with open('rs429358Common.txt', 'r') as rs112:
    #  read the rs7412 site file to concatenate
    with open('rs7412Common.txt', 'r') as rs158:
        # write into a new file
        with open('APE_Merge.txt', 'w') as Merge:
            for line in rs112:
                # removes the newline character at the end of the rs429358 site file
                Merge.write(line.strip('\r\n'))
                # add a separator between the two files
                Merge.write('\t')
                # write the second file after the separator
                Merge.write(rs158.readline())
# identify the APOE genotype
text = open("APE_Merge.txt", "r", encoding='utf-8')
for line in text:
    linesplit = line.split("\t")
    Sample = linesplit[0]
    APE4carrierBefore = linesplit[7]
    C1 = linesplit[8]
    T1 = linesplit[9]
    APE4carrierAfter = linesplit[10]
    C2 = linesplit[18]
    T2 = linesplit[19]
    # assigns a null character to APOEGenetype variables
    APOEGenetype = ""
    # determine the APOE genotype
    if int(C1) == 0 and int(T1) != 0 and int(C2) == 0 and int(T2) != 0:
        APOEGenetype = "E2/2"
    elif int(C1) == 0 and int(T1) != 0 and int(C2) != 0 and int(T2) != 0:
        APOEGenetype = "E2/3"
    elif int(C1) != 0 and int(T1) != 0 and int(C2) != 0 and int(T2) != 0:
        APOEGenetype = "E2/4"
    elif int(C1) == 0 and int(T1) != 0 and int(C2) != 0 and int(T2) == 0:
        APOEGenetype = "E3/3"
    elif int(C1) != 0 and int(T1) != 0 and int(C2) != 0 and int(T2) == 0:
        APOEGenetype = "E3/4"
    elif int(C1) != 0 and int(T1) == 0 and int(C2) != 0 and int(T2) == 0:
        APOEGenetype = "E4/4"
    elif int(C1) == 0 and int(T1) == 0 and int(C2) == 0 and int(T2) == 0:
        APOEGenetype = "NoDetect"
    else:
        APOEGenetype = "CannotJudge"
    output1 = "\t".join(linesplit)
    # remove the newline character
    output1 = output1.replace('\n', '').replace('\r', '')
    # add APOEGenetype column to the original file
    output2 = "\t".join([output1, APOEGenetype])
    # write result into a file
    fp = open('E:/AD_Patient/APE/APE_MergeIdentify.txt', 'a+')
    print(output2, file=fp)
    fp.close()
count = len(open("APE_MergeIdentify.txt", 'r').readlines())

file = open("APE_MergeIdentify.txt", "r", encoding='utf-8')
# check the number of lines
print("-----------Predict APOE Genotype-----------")
n = int(5)
print("Show Result of Prediction:")
for i in range(n):
    print(file.readline().strip())
print("Number Lines of Prediction File:" + str(count))
print("", end='\n')
endtime = time.perf_counter()
print("End of Program Execution")
print("The Total Predicted Time is:", (endtime - start), "s")

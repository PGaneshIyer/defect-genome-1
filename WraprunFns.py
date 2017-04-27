import os, sys
import re
import shutil
import tarfile
import getopt
import warnings

class WraprunFns:
    def walklevel(self, some_dir, level=1):
        """
        ADAPTED FROM: https://stackoverflow.com/questions/229186/os-walk-without-digging-into-directories-below
        LEARNING: Cannot use 'return' with argument inside generator
        """
        some_dir = some_dir.rstrip(os.path.sep)
        assert os.path.isdir(some_dir)
        num_sep = some_dir.count(os.path.sep)
        for root, dirs, files in os.walk(some_dir):
            yield root, dirs, files
            num_sep_this = root.count(os.path.sep)
            if num_sep + level <= num_sep_this:
                del dirs[:]

    def get_subdirs(self, dirPath, minDepth=1, maxDepth=1):
        if (minDepth > 0):
            removeSubDirs = [x[0] for x in self.walklevel(dirPath, minDepth-1)
                             if os.path.isdir(x[0])]
            removeSubDirs = list(set(removeSubDirs))
        else:
            removeSubDirs = []
        allSubDirs = [x[0] for x in self.walklevel(dirPath, maxDepth)
                      if os.path.isdir(x[0])]
        subDirsInRange = [x for x in allSubDirs if x not in removeSubDirs]
        return subDirsInRange

    def remove_items_with_strs(self, dataList, strNeglectList):
        outputList = []
        if (strNeglectList is not None):
            for val in dataList:
                includeFlag = True
                for strNeglect in strNeglectList:
                    if strNeglect in val:
                        includeFlag = False
                        break
                if (includeFlag is True):
                    outputList.append(val)
        else:
            outputList = dataList
        return outputList

    def usage(self):
        # TODO: Implement an actual usage function
        print "Please look at Source Code For Implementation"

    def argParse(self, argv):
        # Adapted from http://www.diveintopython.net/scripts_and_streams/command_line_arguments.html
        argDir = {"input": 'wraprun.json'}
        try:
            opts, args = getopt.getopt(argv, "hd:i:",
                                       ["help", "dir=", "inFile="])
        except getopt.GetoptError:
            self.usage()
            sys.exit(2)
        for opt, arg in opts:
            if opt in ("-h", "--help"):
                self.usage()
                sys.exit()
            elif opt in ("-d", "--dir"):
                argDir["workdir"] = arg
            elif opt in ("-i", "--input"):
                argDir["input"] = arg
        return argDir

    def make_tarfile(self, dirPath, inFileList,
                     tarFileName, option=None):
        with tarfile.open(os.path.join(dirPath, tarFileName),
                          "w:gz") as tar:
            for inFile in inFileList:
                try:
                    tar.add(os.path.join(dirPath, inFile))
                except OSError:
                    warnings.warn("Output File not available")
                    print("FileName: ", inFile)

    def vasp_backup_for_restart(self, dirPath):
        if os.path.isfile(os.path.join(dirPath,
                                       'OUTCAR')):
            tarFileList = [x for x in os.listdir(dirPath)
                           if x.endswith(".tgz")]
            print tarFileList
            if len(tarFileList) > 0:
                restartTmpVal = [re.findall(r'\d+', x)[-1]
                                 for x in tarFileList]
                restartTmpVal.sort(key=float)
                restartVal = int(restartTmpVal[-1])+1
            else:
                restartVal = 0
            tarFileName = 'vaspRunNo_'+str(restartVal)+'.tgz'
            inFileList = ['INCAR', 'POSCAR', 'KPOINTS', 'CONTCAR',
                          'XDATCAR', 'CHGCAR', 'CHG', 'OUTCAR',
                          'DOSCAR', 'PROCAR', 'OSZICAR']
            self.make_tarfile(dirPath, inFileList, tarFileName)
            shutil.copy2(os.path.join(dirPath, 'CONTCAR'),
                         os.path.join(dirPath, 'POSCAR'))

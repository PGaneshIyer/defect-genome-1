import os, sys


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
        for val in dataList:
            includeFlag = True
            for strNeglect in strNeglectList:
                if strNeglect in val:
                    includeFlag = False
                    break
            if (includeFlag is True):
                outputList.append(val)
        return outputList

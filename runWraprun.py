#!/usr/local/bin/python

################################################################
# runWraprun.py creates wraprun tasks that can run in Titan at OLCF.
# This code requires an input json file with the necessary parameters.
# A sample input named "wraprun.json" is provided in the repository.
# This code requires WraprunFns module which is available in defect-genome repository.
# The defect-genome repository is located at: https://github.com/rambalachandran/defect-genome
###############################################################
# USAGE: python runWraprun.py -d <dirPath> -i <jsonFileName>
###############################################################
# Author: Janakiraman Balachandran
################################################################
# TODO: (1) For VASP jobs, check in the dirs if the calculation has converged (with pymatgen).
# TODO: (2) If Converged,
# TODO: (2a) remove the directory from dirList.
# TODO: (2b) Parse Output and store it in database (with pymatgen).
# TODO: (3) If NOT converged,
# TODO: (3a) archive old files,
# TODO: (3b) Create new input files based on type of system, level of theory, no.of steps, force, electronic Convergence
################################################################

import os, sys
import json
import wraprun
from WraprunFns import WraprunFns

if __name__ == '__main__':
    # FIXME: Must readin from a file created by runFireworks.py
    # fw_id_list = [1, 3, 5]
    bundle = wraprun.Wraprun(nocopy=True)
    wraprunFns = WraprunFns()
    argDir = wraprunFns.argParse(sys.argv[1:])
    workdir = argDir["workdir"]
    inFileName = argDir["input"]
    print("RootDir: ", workdir)
    print("Input File Name: ", inFileName)
    with open(os.path.join(workdir, inFileName)) as inFile:
        data = json.load(inFile)
    # PRINT FOR VALIDATION
    # for key, value in data.iteritems():
    #     print key, value
    if (data['DirList']):
        dirList = data['DirList']
    else:
        dirListFull = wraprunFns.get_subdirs(workdir,
                                             data['MinDirDepth'],
                                             data['MaxDirDepth'])
        dirList = wraprunFns.remove_items_with_strs(dirListFull,
                                                    data["StringNeglectList"])
    print "VaspData Backup: ", data["VaspBackup"]
    if (data["VaspBackup"] == 1):
        for dirVal in dirList:
            print dirVal
            wraprunFns.vasp_backup_for_restart(dirVal)
    else:
        if (data["ProcType"] == "CPU"):
            if (data["ProcsPerTaskList"]):
                taskList = data["ProcsPerTaskList"]
            else:
                taskList = [data['ProcsPerTask']]*len(dirList)
            pes_per_node = data['ProcsPerNode']
            bundle.add_task(pes=taskList, pes_per_node=pes_per_node,
                            exe=data['Executable'], cd=dirList)
        elif (data["ProcType"] == "GPU"):
            if (data["NodesPerTaskList"]):
                taskList = data["NodesPerTaskList"]
            else:
                taskList = [data['NodesPerTask']]*len(dirList)
            pes_per_node = data['TasksPerNode']
            bundle.add_task(pes=taskList, pes_per_node=pes_per_node,
                            depth=data["ThreadsPerTask"],
                            exe=data['Executable'], cd=dirList)
        print("No. of Dirs: ", len(dirList))
        print("DirList: ", dirList)
        bundle.launch()

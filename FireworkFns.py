import os,sys
from fireworks import Firework, LaunchPad, ScriptTask, Workflow, Launch
from fireworks.core.rocket_launcher import launch_rocket
from fireworks_vasp.tasks import WriteVaspInputTask, VaspCustodianTask, VaspAnalyzeTask
from pymatgen.core.structure import Structure
from pymatgen.io.vaspio import Poscar

launchpad = LaunchPad(host='localhost', port=2221, name='*',
                      username='*', password='*')

def get_fw_info(fw_id, specDictKeyList):
    specDict = {}
    firework = launchpad.get_fw_by_id(fw_id)
    for key in specDictKeyList:
        specDict[key] = firework.spec[key]
    return firework.name, firework.state, firework.launches[-1].launch_id, firework.launches[-1].launch_dir, specDict


def get_crystal_struc(launch_id, getInputStruc=False):
    launch = launchpad.get_launch_by_id(launch_id)
    if getInputStruc:
        if launch.action:
            strucDict = launch.action.stored_data.get('vasprun', {}).get('input',{}).get('crystal', {})
            if len(strucDict) > 0:
                crystalStruc = Structure.from_dict(strucDict)
            else:
                strucFilePath = os.path.join(launch.launch_dir, 'POSCAR')
                crystalStruc = Structure.from_file(strucFilePath)
        else:
            strucFilePath = os.path.join(launch.launch_dir, 'POSCAR')
            crystalStruc = Structure.from_file(strucFilePath)

    else:
        if launch.action:
            strucDict = launch.action.stored_data.get('vasprun',
                                                       {}).get('output', {}).get('crystal', {})
            if len(strucDict) > 0:
                crystalStruc = Structure.from_dict(strucDict)
            else:
                strucFilePath = os.path.join(launch.launch_dir, 'CONTCAR')
                crystalStruc = Structure.from_file(strucFilePath)
        else:
            strucFilePath = os.path.join(launch.launch_dir, 'CONTCAR')
            crystalStruc = Structure.from_file(strucFilePath)
    return crystalStruc

def create_fireworks(structure, name, specDict, fw_id,
                     viset='CubicGGA24PVInputSet',                params={"user_incar_settings":{"NPAR":"16","NSIM":"4","LXML":"TRUE"}}, handlers=["DetourErrorHandler","PBSWalltimeHandler"], vasp_cmd=["aprun","-n","128","vasp"]):
    name = name
    t1 = WriteVaspInputTask(structure=structure, vasp_input_set=viset, input_set_params=params)
    t2 = VaspCustodianTask(vasp_cmd=vasp_cmd, handlers=handlers)
    t3 = VaspAnalyzeTask()
    firework = Firework([t1, t2, t3], name=name, spec=specDict)
    return firework

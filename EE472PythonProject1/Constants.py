from enum import Enum

class StateOfReader(Enum):
    EmpthyReadingMode = 0
    BusReadingMode = 1
    BranchReadingMode = 2

class TypeOfBara(Enum):
    Slack = 0
    PV = 1
    PQ = 2
class DataOfBusForSimulation:
    def __init__(self,
                 busnumber: int,
                 busindex: int,
                 type: TypeOfBara,
                 finalvoltage: float,
                 finalangle: float,
                 loadp: float,
                 loadq: float,
                 generationp: float,
                 generationq: float,
                 qorvmaxlimit: float,
                 qorvminlimit: float):
        self.BusIndex = busindex
        self.BusNumber = busnumber
        self.Type = type
        self.FinalVoltage = finalvoltage
        self.FinalAngle = finalangle
        self.LoadP = loadp
        self.LoadQ = loadq
        self.GenerationP = generationp
        self.GenerationQ = generationq
        self.QorVMaxLimit = qorvmaxlimit
        self.QorVMinLimit = qorvminlimit


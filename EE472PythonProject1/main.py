import DataReaderandYBUSCreater as R
import NRSimulation as S
import ControlSimulation as C
import time

import numpy as np
import math as mat


def Run(directory,MAXITERATION, EPSILON):
    # reading data and creating YBUS also the type of the baras
    (YBUS, SlackBara, PVBaras, PQBaras) = R.Read(directory)
    # running the Newton Raphsen method
    (currentVoltages, currentAngles) = S.Simulate(YBUS, SlackBara, PVBaras, PQBaras, MAXITERATION, EPSILON, False)
    # controlling the datas obtained from the simulation
    C.Control(currentVoltages, currentAngles, SlackBara, PVBaras, PQBaras)

    TotalP = 0
    TotalQ = 0

    for PVbara in PVBaras:
        TotalP = TotalP + PVbara.GenerationP - PVbara.LoadP
        TotalQ = TotalQ + PVbara.GenerationQ - PVbara.LoadQ
    for PQBara in PQBaras:
        TotalP = TotalP + PQBara.GenerationQ - PQBara.LoadP
        TotalQ = TotalQ + PQBara.GenerationQ - PQBara.LoadQ

    currentIndexBus = SlackBara.BusIndex

    # taking the current voltage of the bus
    currentVoltage = currentVoltages[currentIndexBus]

    # initializing the current P and Q
    Pi = 0
    Qi = 0

    for busIndexColumn in range(118):
        # taking the data from YBus
        Yik = YBUS[currentIndexBus, busIndexColumn]

        # reel part of the YBUS
        Gik = np.real(Yik)

        # imaginery part of YBUS
        Bik = np.imag(Yik)

        # angle difference
        thetaik = currentAngles[currentIndexBus] - currentAngles[busIndexColumn]

        Pi = Pi + currentVoltage * currentVoltages[busIndexColumn] * (
                Gik * mat.cos(thetaik[0]).real + Bik * mat.sin(thetaik[0]).real)
        Qi = Qi + currentVoltage * currentVoltages[busIndexColumn] * (
                Gik * mat.sin(thetaik[0]).real - Bik * mat.cos(thetaik[0]).real)


    print("Slack bara Pi calculation: " + str(Pi + SlackBara.LoadP) + "Slack bara Qi calculation: " + str(Qi))
    print("-----------------------------------------------------")
    print("TotalP: " + str(TotalP + float(Pi + SlackBara.LoadP)) + " TotalQ: " + str(TotalQ + float(Qi)))
    print("-----------------------------------------------------")


if __name__ == "__main__":
    start_time = time.time()

    Run("ieee118cdf.txt", 25, 5e-13)

    end_time = time.time()

    execution_time = end_time - start_time

    print("Execution time:", execution_time, "seconds")


















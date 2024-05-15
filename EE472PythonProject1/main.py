
import time
import matplotlib.pyplot as plt
import cmath as mat
import numpy as np
import math as mat
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



def Print(currentVoltages, currentAngles, solution_time, numberOfIterations, realLoss, reactiveLoss, PVBaras, PQBaras, SlackBara):

    print("-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------")
    print("------------------------------------------------------------------------------------------------OUTPUT DATA------------------------------------------------------------------------------------------------")
    print("-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------")

    for index in range(len(currentVoltages)):

        isBaraDataFound = False

        if(SlackBara.BusIndex == index):
            baraData = SlackBara
            isBaraDataFound = True

        if(not isBaraDataFound):
            for Data in PVBaras:
                if (Data.BusIndex == index):
                    baraData = Data
                    isBaraDataFound = True
                    break

        if (not isBaraDataFound):
            for Data in PQBaras:
                if (Data.BusIndex == index):
                    baraData = Data
                    isBaraDataFound = True
                    break

        print("Bus Number: " + str(baraData.BusNumber) + " Voltages: (in pu): " + str(currentVoltages[index][0]) + " Voltage angles (in degree): " + str(currentAngles[index][0]*180/mat.pi) + " ----> " + str(baraData.Type))

    print("-----------------------------------------------")
    print("Solution Time in (sec): " + str(solution_time))
    print("Number of iterations: " + str(numberOfIterations))
    print("Real Power Loss (in pu): " + str(realLoss))
    print("Reactive Power Loss (in pu): " + str(reactiveLoss))
    print(
        "-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------")
    print(
        "-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------")
    print(
        "-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------")


def CalculateLoss(SlackBara ,PVBaras, PQBaras, currentVoltages, currentAngles, YBUS):
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


    for busIndexColumn in range(len(PVBaras)+len(PQBaras)+1):
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


    # print("Total Real Power Loss: " + str(TotalP + float(Pi + SlackBara.LoadP)))
    # print("Total Reactive Power Loss: " + str(TotalQ + float(Qi + SlackBara.LoadQ)))

    return TotalP + float(Pi + SlackBara.LoadP), TotalQ + float(Qi + SlackBara.LoadQ)
def Control(currentVoltages, currentAngles, SlackBara, PVBaras, PQBaras):
    ### creating control mechanism for the simulation ###
    busCounter = len(PVBaras) + len(PQBaras) + 1

    TrueVoltages = np.zeros((busCounter, 1))
    TrueAngles = np.zeros((busCounter, 1))

    TrueVoltages[SlackBara.BusIndex] = SlackBara.FinalVoltage
    TrueAngles[SlackBara.BusIndex] = SlackBara.FinalAngle

    for PVbara in PVBaras:
        TrueVoltages[PVbara.BusIndex] = PVbara.FinalVoltage
        TrueAngles[PVbara.BusIndex] = PVbara.FinalAngle
    for PQbara in PQBaras:
        TrueVoltages[PQbara.BusIndex] = PQbara.FinalVoltage
        TrueAngles[PQbara.BusIndex] = PQbara.FinalAngle

    # flag
    IsThereError = False

    # Compute the absolute difference between the column vectors
    differenceVoltages = np.abs(currentVoltages - TrueVoltages)

    # Find indices where the difference exceeds 1e-10
    indices = np.where(differenceVoltages > 1e-2)

    # Print the indices where the difference exceeds the threshold voltages
    print("Indices where the difference exceeds (for voltages) 1e-2:")
    for index in zip(indices[0]):
        IsThereError = True
        print("index: " + str(index) + " calculated voltage: " + str(currentVoltages[index]) + " real voltage: " + str(
            TrueVoltages[index]))

    # Compute the absolute difference between the column vectors
    differenceAngles = np.abs(currentAngles - TrueAngles) * 180 / mat.pi

    # Find indices where the difference exceeds 1e-10
    indices = np.where(differenceAngles > 1e-0)

    # Print the indices where the difference exceeds the threshold voltages
    print("Indices where the difference exceeds (for angles) 1e-0:")
    for index in zip(indices[0]):
        IsThereError = True
        print("index: " + str(index) + " calculated angle: " + str(
            currentAngles[index] * 180 / mat.pi) + " real angle: " + str(TrueAngles[index] * 180 / mat.pi))

    if (not IsThereError):
        print("SUCCESS")
    else:
        print("FAILURE")
    #####################################################
    # print(currentVoltages)
    # print(currentAngles*180/mat.pi)
    # print(TrueAngles*180/mat.pi)

def Simulate(YBUS, SlackBara, PVBaras, PQBaras, MAXITERATION, EPSILON, HotStartFlag):
    ################################################################################ Newton Raphsen method ################################################################################
    # flat start conditions
    flatStartVoltageMagnitude = SlackBara.FinalVoltage
    flatStartAngle = SlackBara.FinalAngle

    busCounter = len(PVBaras) + len(PQBaras) + 1
    epsilon = EPSILON

    # creating the voltage and angle vectors
    currentVoltages = np.ones((busCounter, 1)) * flatStartVoltageMagnitude
    currentAngles = np.ones((busCounter, 1)) * flatStartAngle

    PVBussesThatAreStuckInQlimits = []

    # assign PV bara voltages from data
    for PVbarainfo in PVBaras:
        currentVoltages[PVbarainfo.BusIndex] = PVbarainfo.FinalVoltage

    Big_Loop_Index = 0

    # enable or disable correct values to angles and voltages
    hotStartBool = HotStartFlag
    if (hotStartBool):
        # setting slack bara data to voltage and angle
        SlackBaraIndex = SlackBara.BusIndex
        SlackBaraVoltage = SlackBara.FinalVoltage
        SlackBaraAngle = SlackBara.FinalAngle
        currentVoltages[SlackBaraIndex] = SlackBaraVoltage
        currentAngles[SlackBaraIndex] = SlackBaraAngle
        # setting PVbaras data to voltage and angle
        for PVbara in PVBaras:
            PVbaraindex = PVbara.BusIndex
            PVbaraVoltage = PVbara.FinalVoltage
            PVbaraAngle = PVbara.FinalAngle
            currentVoltages[PVbaraindex] = PVbaraVoltage
            currentAngles[PVbaraindex] = PVbaraAngle
        # setting PQbaras data to voltage and angle
        for PQbara in PQBaras:
            PQbaraindex = PQbara.BusIndex
            PQbaraVoltage = PQbara.FinalVoltage
            PQbaraAngle = PQbara.FinalAngle
            currentVoltages[PQbaraindex] = PQbaraVoltage
            currentAngles[PQbaraindex] = PQbaraAngle

    MAXITERATIONExceededInfoFlag = False
    while (Big_Loop_Index < MAXITERATION):
        print("Index: " + str(Big_Loop_Index))
        Big_Loop_Index = Big_Loop_Index + 1
        if(Big_Loop_Index == MAXITERATION):
            MAXITERATIONExceededInfoFlag = True
        ################ Checking constraints, calculating P and Q values of all current baras and making transfers between PV to PQ baras ################
        # initializing Pi and Qi data containers
        # For PV and PQ baras both Pi and Qi values are calculated
        Pis = np.zeros((len(PVBaras) + len(PQBaras), 1))
        Qis = np.zeros((len(PVBaras) + len(PQBaras), 1))

        # calculation for PVbaras P and Q
        index = 0
        listOfbarasPVtoPQ = []
        for dataForBus in PVBaras:

            # taking current bus index
            currentIndexBus = dataForBus.BusIndex

            # taking the current voltage of the bus
            currentVoltage = currentVoltages[currentIndexBus]

            # initializing the current P and Q
            Qi = 0
            Pi = 0

            for busIndexColumn in range(busCounter):
                # taking the data from YBus
                Yik = YBUS[currentIndexBus, busIndexColumn]

                # reel part of the YBUS
                Gik = np.real(Yik)

                # imaginery part of YBUS
                Bik = np.imag(Yik)

                # angle difference
                thetaik = currentAngles[currentIndexBus] - currentAngles[busIndexColumn]

                Pi = Pi + currentVoltage * currentVoltages[busIndexColumn] * (
                            Gik * mat.cos(thetaik[0]).real + Bik * mat.sin(thetaik[0])).real
                Qi = Qi + currentVoltage * currentVoltages[busIndexColumn] * (
                            Gik * mat.sin(thetaik[0]).real - Bik * mat.cos(thetaik[0])).real

            # checking reactive power limits
            # all the baras in PV baras has limits
            # if(dataForBus.QorVMaxLimit < Qi + dataForBus.LoadQ or dataForBus.QorVMinLimit > Qi + dataForBus.LoadQ):
            #     # add the PV bara to the transfer list if its Q value is on the outside of boundries
            #     listOfbarasPVtoPQ.append(dataForBus)
            #     PVBussesThatAreStuckInQlimits.append(dataForBus)
            #     # since bara will be transferred index will not be incremented (no data stored)
            #     continue
            # else:
            #     # if PV bara satisfy Q boundries store its P and Q data
            #     Pis[index] = Pi
            #     Qis[index] = Qi
            #     index = index + 1

            Pis[index] = Pi
            Qis[index] = Qi
            index = index + 1
            # print("bus index: " + str(dataForBus.BusIndex))
            # print("MAXQ: " + str(dataForBus.QorVMaxLimit))
            # print("MINQ: " + str(dataForBus.QorVMinLimit))
            # print("currentQ: " + str(Qi))
            # print("currentV: " + str(currentVoltages[dataForBus.BusIndex]))
            # print("finalV: " + str(dataForBus.FinalVoltage))
            # Qis[index] = Qi
            # Pis[index] = Pi
            # index = index + 1

        # if the Q calculation remain outside the constraints make PV bara PQ
        # for bara in listOfbarasPVtoPQ:
        #     PQBaras.append(bara)
        #     PVBaras.remove(bara)



        # calculation for PQbaras P and Q
        for dataForBus in PQBaras:

            # taking current bus index
            currentIndexBus = dataForBus.BusIndex

            # taking the current voltage of the bus
            currentVoltage = currentVoltages[currentIndexBus]

            # initializing the current P and Q
            Pi = 0
            Qi = 0

            for busIndexColumn in range(busCounter):
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

            # checking reactive power limits
            if(dataForBus.Type == TypeOfBara.PV):
                if(dataForBus.QorVMaxLimit < Qi + dataForBus.LoadQ):
                    # if PQ baras Q is at limits saturate the value (max limit)
                    # print("bus index: " + str(dataForBus.BusIndex))
                    # print("MAXQ: " + str(dataForBus.QorVMaxLimit))
                    # print("currentQ: " + str(Qi))
                    # print("currentV: " + str(currentVoltages[dataForBus.BusIndex]))
                    # print("finalV: " + str(dataForBus.FinalVoltage))

                    # print("before saturation Qi (max) (bus id: " + str(dataForBus.BusIndex + 1) + "): " + str(Qi))
                    # print("before saturation Qig (max) (bus id: " + str(dataForBus.BusIndex + 1) + "): " + str(Qi + dataForBus.LoadQ))
                    # Qi = dataForBus.QorVMaxLimit - dataForBus.LoadQ
                    # print("after saturation Qi (max) (bus id: " + str(dataForBus.BusIndex + 1) + "): " + str(Qi))
                    # print("after saturation Qig (max) (bus id: " + str(dataForBus.BusIndex + 1) + "): " + str(Qi + dataForBus.LoadQ))
                    pass
                elif (dataForBus.QorVMinLimit > Qi + dataForBus.LoadQ):
                    # if PQ baras Q is at limits saturate the value (min limit)
                    # print("bus index: " + str(dataForBus.BusIndex))
                    # print("MAXQ: " + str(dataForBus.QorVMinLimit))
                    # print("currentQ: " + str(Qi))
                    # print("currentV: " + str(currentVoltages[dataForBus.BusIndex]))
                    # print("finalV: " + str(dataForBus.FinalVoltage))

                    # print("before saturation Qi (min) (bus id: " + str(dataForBus.BusIndex+1) + "): " + str(Qi))
                    # print("before saturation Qig (min) (bus id: " + str(dataForBus.BusIndex + 1) + "): " + str(Qi + dataForBus.LoadQ))
                    # Qi = dataForBus.QorVMinLimit - dataForBus.LoadQ
                    # print("after saturation Qi (min) (bus id: " + str(dataForBus.BusIndex + 1) + "): " + str(Qi))
                    # print("after saturation Qig (min) (bus id: " + str(dataForBus.BusIndex + 1) + "): " + str(Qi + dataForBus.LoadQ))
                    pass



            Pis[index] = Pi
            Qis[index] = Qi

            index = index + 1

        # for bara in PVBussesThatAreStuckInQlimits:
        #     indexInPQBara = PQBaras.index(bara)
        #     indexInQis = indexInPQBara + len(PVBaras)
        #     if(abs(Qis[indexInQis] + bara.LoadQ - bara.QorVMaxLimit) < 1e-8):
        #         print("Bus Id: " + str(bara.BusIndex + 1) + " stuck at max limit")
        #     elif(abs(Qis[indexInQis] + bara.LoadQ - bara.QorVMinLimit) < 1e-8):
        #         print("Bus Id: " + str(bara.BusIndex + 1) + " stuck at min limit")

        # if actual PV bara in PQ baras return back to its limits make it PV bara
        # listOfBarasPQtoPV = []
        # listOfIndexOfTransferredBaras = []
        # index = 0
        # for dataForBus in PQBaras:
        #     if (dataForBus.Type == const.TypeOfBara.PV):
        #         # for busIndexColumn in range(busCounter - 1):
        #         #     # taking the data from YBus
        #         #     Yik = selfCreatedYBusMatrix57[currentIndexBus, busIndexColumn]
        #         #     # reel part of the YBUS
        #         #     Gik = np.real(Yik)
        #         #     # imaginery part of YBUS
        #         #     Bik = np.imag(Yik)
        #         #     # angle difference
        #         #     thetaik = currentAngles[currentIndexBus] - currentAngles[busIndexColumn]
        #         #     Qi = Qi + currentVoltage * currentVoltages[busIndexColumn] * (Gik * mat.sin(thetaik[0]).real - Bik * mat.cos(thetaik[0]).real)
        #         currentIndexBus = dataForBus.BusIndex
        #         currentVoltage = currentVoltages[currentIndexBus]
        #         QiOfPVBaraInPQbara = Qis[index + len(PVBaras)] + dataForBus.LoadQ
        #         # print("bus index:")
        #         # print(dataForBus.BusIndex)
        #         # print("current Qi:")
        #         # print(QiOfPVBaraInPQbara)
        #         # print("current voltage")
        #         # print(currentVoltage)
        #         # print("final voltage")
        #         # print(dataForBus.FinalVoltage)
        #         if (QiOfPVBaraInPQbara > dataForBus.QorVMinLimit and QiOfPVBaraInPQbara < dataForBus.QorVMaxLimit):
        #             listOfBarasPQtoPV.append(dataForBus)
        #             listOfIndexOfTransferredBaras.append(index)
        #     index = index + 1
        # # adjusting positions of the Pis and Qis since we are changing positions of the PV and PQ baras
        # PisTemp = np.zeros((len(PVBaras) + len(PQBaras), 1))
        # QisTemp = np.zeros((len(PVBaras) + len(PQBaras), 1))
        # # print("length of pv baras")
        # # print(len(PVBaras))
        # # print("length of pq baras")
        # # print(len(PQBaras))
        # # print("length of transferred baras")
        # # print(len(listOfIndexOfTransferredBaras))
        # indeces = []
        # lastPartIndex = 0
        # for index in range(len(PVBaras) + len(PQBaras)):
        #     # first add PV bara datas
        #     if(index < len(PVBaras)):
        #         PisTemp[index] = Pis[index]
        #         QisTemp[index] = Qis[index]
        #     # second add appended PV bara datas
        #     elif(index < len(PVBaras) + len(listOfIndexOfTransferredBaras)):
        #         littleIndex = index - len(PVBaras)
        #         currentPowerIndex = len(PVBaras) + listOfIndexOfTransferredBaras[littleIndex]
        #         indeces.append(currentPowerIndex)
        #         PisTemp[index] = Pis[currentPowerIndex]
        #         QisTemp[index] = Qis[currentPowerIndex]
        #     #third add remaining data
        #     else:
        #         while(indeces.__contains__(index - len(listOfIndexOfTransferredBaras) + lastPartIndex)):
        #             lastPartIndex = lastPartIndex + 1
        #         littleIndex = index - len(listOfIndexOfTransferredBaras) + lastPartIndex
        #         PisTemp[index] = Pis[littleIndex]
        #         QisTemp[index] = Qis[littleIndex]
        # # print(np.concatenate((QisTemp,Qis),axis=1))
        # Pis = PisTemp
        # Qis = QisTemp
        # # if the voltage constraint violated sent PQ bara to PV again (old PV bara)
        # for bara in listOfBarasPQtoPV:
        #     PQBaras.remove(bara)
        #     PVBaras.append(bara)
        #     # currentVoltages[bara.BusIndex] = bara.FinalVoltage
        #     PVBussesThatAreStuckInQlimits.remove(bara)



        ###################################################################################################################################################

        ################ calculating the deltaP for PV and PQ baras and deltaQ for PQ baras ################
        # initializing the deltaP and deltaQ vector
        # for deltaP first PV Baras then PQ baras are stored
        # for deltaQ only the PQ baras are stored
        deltaPs = np.zeros((len(PVBaras) + len(PQBaras), 1))
        DeltaQs = np.zeros((len(PQBaras), 1))
        Deltas = np.concatenate((deltaPs, DeltaQs), axis=0)

        # index of delta vector and Pi vector is same (index -> index of the delta vector)
        index = 0
        # calculating the deltaPi for PV Baras
        for dataForBus in PVBaras:
            # calculating the target P
            Pispec = dataForBus.GenerationP - dataForBus.LoadP

            deltaPi = Pispec - Pis[index]
            Deltas[index] = deltaPi
            index = index + 1

        # calculating the deltaPi for PQ Baras
        for dataForBus in PQBaras:
            # calculating the target P
            Pispec = dataForBus.GenerationP - dataForBus.LoadP

            deltaPi = Pispec - Pis[index]
            Deltas[index] = deltaPi
            index = index + 1

        # index of delta vector and Qi vector is not same
        indexQiVector = len(PVBaras)
        # calculating the deltaQi
        for dataForBus in PQBaras:
            # calculating the target Q
            Qispec = dataForBus.GenerationQ - dataForBus.LoadQ

            deltaQi = Qispec - Qis[indexQiVector]
            Deltas[index] = deltaQi
            index = index + 1
            indexQiVector = indexQiVector + 1
        ####################################################################################################

        ############################# STOPPING CONDITION OF THE SIMULATION #############################
        maxDif = np.max(np.abs(Deltas))
        print("Max difference deltaP and deltaQ: " + str(maxDif))
        if (maxDif < epsilon):
            break
        ################################################################################################

        ########################################################## CONSTRUCTING BIG MATRIX ##########################################################
        tempPVPQConcatenated = PVBaras + PQBaras

        ################ constructing the H matrix ################
        size_of_H_matrix = len(PVBaras) + len(PQBaras)
        H = np.zeros((size_of_H_matrix, size_of_H_matrix))

        for y in range(size_of_H_matrix):

            i = tempPVPQConcatenated[y].BusIndex

            for x in range(size_of_H_matrix):

                j = tempPVPQConcatenated[x].BusIndex

                # if i == j then x == y for H matrix
                if (i == j):
                    # finding Qi
                    Qi = Qis[y]

                    # taking the data from YBus
                    Yii = YBUS[i, i]

                    # imaginery part of YBUS
                    Bii = np.imag(Yii)

                    # current voltage
                    Vi = currentVoltages[i]

                    H[y, x] = (-1) * Qi - Bii * Vi * Vi
                elif (i != j):
                    # finding thetaij
                    thetaij = currentAngles[i] - currentAngles[j]

                    # finding Vi
                    Vi = currentVoltages[i]

                    # finding Vj
                    Vj = currentVoltages[j]

                    # taking the data from YBus
                    Yij = YBUS[i, j]

                    # imaginery part of YBUS
                    Bij = np.imag(Yij)

                    # real part of YBUS
                    Gij = np.real(Yij)

                    H[y, x] = Vi * Vj * (Gij * mat.sin(thetaij[0]).real - Bij * mat.cos(thetaij[0]).real)
        ###########################################################

        ################ constructing the N matrix ################
        number_of_columns_N = len(PQBaras)
        number_of_rows_N = len(PQBaras) + len(PVBaras)
        N = np.zeros((number_of_rows_N, number_of_columns_N))

        for y in range(number_of_rows_N):

            i = tempPVPQConcatenated[y].BusIndex

            for x in range(number_of_columns_N):

                j = PQBaras[x].BusIndex

                # if i == j
                if (i == j):
                    # finding Pi
                    Pi = Pis[y]

                    # taking the data from YBus
                    Yii = YBUS[i, i]

                    # real part of YBUS
                    Gii = np.real(Yii)

                    # current voltage
                    Vi = currentVoltages[i]

                    N[y, x] = Pi + Gii * Vi * Vi
                elif (i != j):
                    # finding thetaij
                    thetaij = currentAngles[i] - currentAngles[j]

                    # finding Vi
                    Vi = currentVoltages[i]

                    # finding Vj
                    Vj = currentVoltages[j]

                    # taking the data from YBus
                    Yij = YBUS[i, j]

                    # imaginery part of YBUS
                    Bij = np.imag(Yij)

                    # real part of YBUS
                    Gij = np.real(Yij)

                    N[y, x] = Vi * Vj * (Gij * mat.cos(thetaij[0]).real + Bij * mat.sin(thetaij[0]).real)
        ###########################################################

        ################ constructing the M matrix ################
        number_of_columns_M = len(PQBaras) + len(PVBaras)
        number_of_rows_M = len(PQBaras)
        M = np.zeros((number_of_rows_M, number_of_columns_M))

        for y in range(number_of_rows_M):

            i = PQBaras[y].BusIndex

            for x in range(number_of_columns_M):

                j = tempPVPQConcatenated[x].BusIndex

                # if i == j
                if (i == j):
                    # finding Pi
                    Pi = Pis[x]

                    # taking the data from YBus
                    Yii = YBUS[i, i]

                    # real part of YBUS
                    Gii = np.real(Yii)

                    # current voltage
                    Vi = currentVoltages[i]

                    M[y, x] = Pi - Gii * Vi * Vi
                elif (i != j):
                    # finding thetaij
                    thetaij = currentAngles[i] - currentAngles[j]

                    # finding Vi
                    Vi = currentVoltages[i]

                    # finding Vj
                    Vj = currentVoltages[j]

                    # taking the data from YBus
                    Yij = YBUS[i, j]

                    # imaginery part of YBUS
                    Bij = np.imag(Yij)

                    # real part of YBUS
                    Gij = np.real(Yij)

                    M[y, x] = (-1) * Vi * Vj * (Gij * mat.cos(thetaij[0]).real + Bij * mat.sin(thetaij[0]).real)
        ###########################################################

        ################ constructing the L matrix ################
        size_of_L_matrix = len(PQBaras)
        L = np.zeros((size_of_L_matrix, size_of_L_matrix))

        for y in range(size_of_L_matrix):

            i = PQBaras[y].BusIndex

            for x in range(size_of_L_matrix):

                j = PQBaras[x].BusIndex

                # if i == j
                if (i == j):
                    # finding Qi
                    Qi = Qis[len(PVBaras) + x]

                    # taking the data from YBus
                    Yii = YBUS[i, i]

                    # imag part of YBUS
                    Bii = np.imag(Yii)

                    # current voltage
                    Vi = currentVoltages[i]

                    L[y, x] = Qi - Bii * Vi * Vi
                elif (i != j):
                    # finding thetaij
                    thetaij = currentAngles[i] - currentAngles[j]

                    # finding Vi
                    Vi = currentVoltages[i]

                    # finding Vj
                    Vj = currentVoltages[j]

                    # taking the data from YBus
                    Yij = YBUS[i, j]

                    # imaginery part of YBUS
                    Bij = np.imag(Yij)

                    # real part of YBUS
                    Gij = np.real(Yij)

                    L[y, x] = Vi * Vj * (Gij * mat.sin(thetaij[0]).real - Bij * mat.cos(thetaij[0]).real)
        ###########################################################

        ################ concatenating the matrices and finalizing the Big_Matrix_Inverse ################
        Left = np.concatenate((H, M), axis=0)
        Right = np.concatenate((N, L), axis=0)

        Big_Matrix = np.concatenate((Left, Right), axis=1)
        Big_Matrix_Inverse = np.linalg.inv(Big_Matrix)
        ##################################################################################################
        #############################################################################################################################################

        ############ finding the theta and voltage differences ############
        DeltasThetaAndVoltage = np.dot(Big_Matrix_Inverse, Deltas)

        first_array = DeltasThetaAndVoltage[:len(PVBaras) + len(PQBaras)]
        second_array = DeltasThetaAndVoltage[len(PVBaras) + len(PQBaras):]

        #print("Max difference angle: " + str(np.max(first_array)) + " Max voltage difference ratio: " + str(np.max(second_array)))
        ###################################################################
        ####### finding the next iteration voltage and theta values #######
        for index in range(DeltasThetaAndVoltage.shape[0]):
            if (index < len(tempPVPQConcatenated)):
                # for delta thetas
                busindex = tempPVPQConcatenated[index].BusIndex
                currentAngles[busindex] = currentAngles[busindex] + DeltasThetaAndVoltage[index]
            else:
                # for delta voltages
                indexInPQBaras = index - len(tempPVPQConcatenated)
                busindex = PQBaras[indexInPQBaras].BusIndex
                currentVoltages[busindex] = currentVoltages[busindex] * (1 + DeltasThetaAndVoltage[index])
        ###################################################################
        # print(Big_Matrix)
        # file_path = "C:/Users/DELL/Desktop/EE472ProjectCorrectYbus/Big_matrix" + str(Big_Loop_Index) + ".mat"
        # savemat(file_path,{'matrix'+ str(Big_Loop_Index):Big_Matrix_Inverse})
        # print(np.linalg.det(Big_Matrix))


    numberOfIteration = Big_Loop_Index + 1

    if(MAXITERATIONExceededInfoFlag):
        print("MAXITERATION EXCEEDED!!!!!!!!!")
    return currentVoltages, currentAngles, numberOfIteration
    #######################################################################################################################################################################################

def Read(cdfDirectory):

    ieeecdf = open(cdfDirectory,"r")

    ########## creating YBus matrix ##########
    readingMode = StateOfReader.EmpthyReadingMode
    busNumberToDataOfBus = {}
    PQBaras = []
    PVBaras = []
    busCounter = 0
    TotalNumberOfBus = 0
    ComplexPowerBase = 0

    previousLineDataHolderForGettingPowerBase = 0


    ########## Reading data and Creating YBUS ##########
    for line in ieeecdf:

        ##### reading mode selection #####
        if (line.__contains__("BRANCH DATA FOLLOWS ")):
            readingMode = StateOfReader.BranchReadingMode
            continue
        elif (line.__contains__("BUS DATA FOLLOWS ")):
            readingMode = StateOfReader.BusReadingMode

            ComplexPowerBase = float(previousLineDataHolderForGettingPowerBase[31:37].strip())

            # Split the string by whitespace
            words = line.split()

            # Iterate through each word to find the numeric member
            for word in words:
                if word.isdigit():  # Check if the word is numeric
                    TotalNumberOfBus = int(word)
                    break

            YBUS = np.zeros((TotalNumberOfBus, TotalNumberOfBus), dtype=np.complex128)

            continue
        elif (line.__contains__("-999")):
            readingMode = StateOfReader.EmpthyReadingMode
            continue

        ##### reading section #####
        ### reading bus data ###
        if (readingMode == StateOfReader.BusReadingMode):
            # debugging the data
            # print("----------------------------------------")
            # print("--bus number--")
            # print(line[0:4].strip())
            bus_number = int(line[0:4].strip())
            # print(bus_number)
            # print("--shunt susceptance--")
            # print(line[114:122].strip())
            bus_susceptance = float(line[114:122].strip())
            # print(bus_susceptance)
            # print("--shunt conductance--")
            # print(line[106:114].strip())
            bus_conductance = float(line[106:114].strip())
            # print(bus_conductance)
            # print("--type of bara--")
            # print(line[24:26].strip())
            type_of_bara_int = int(line[24:26].strip())
            # print(type_of_bara_int)
            type_of_bara = TypeOfBara.Slack
            if (type_of_bara_int == 0 or type_of_bara_int == 1):
                type_of_bara = TypeOfBara.PQ
            elif (type_of_bara_int == 2):
                type_of_bara = TypeOfBara.PV



            # print("--final voltage--")
            # print(line[27:33].strip())
            final_voltage = float(line[27:33].strip())
            # print(final_voltage)
            # print("--final angle--")
            # print(line[33:40].strip())
            final_angle = float(line[33:40].strip()) * mat.pi / 180  # in radians
            # print(final_angle)
            # print("--load MW--")
            # print(line[40:49].strip())
            load_MW = float(line[40:49].strip()) / ComplexPowerBase
            # print(load_MW)
            # print("--load MVAR--")
            # print(line[49:59].strip())
            load_MVAR = float(line[49:59].strip()) / ComplexPowerBase
            # print(load_MVAR)
            # print("--generation MW--")
            # print(line[59:67].strip())
            generation_MW = float(line[59:67].strip()) / ComplexPowerBase
            # print(generation_MW)
            # print("--generation MVAR--")
            # print(line[67:75].strip())
            generation_MVAR = float(line[67:75].strip()) / ComplexPowerBase
            # print(generation_MVAR)
            # print("--max MVAR or voltage limit--")
            # print(line[90:98].strip())
            max_MVAR_or_voltage_limit = float(line[90:98].strip()) / ComplexPowerBase
            # print(max_MVAR_or_voltage_limit)
            # print("--min MVAR or voltage limit--")
            # print(line[98:106].strip())
            min_MVAR_or_voltage_limit = float(line[98:106].strip()) / ComplexPowerBase
            # print(min_MVAR_or_voltage_limit)
            # print("----------------------------------------")

            current_bus_data = DataOfBusForSimulation(bus_number, busCounter, type_of_bara, final_voltage,
                                                            final_angle, load_MW, load_MVAR, generation_MW,
                                                            generation_MVAR, max_MVAR_or_voltage_limit,
                                                            min_MVAR_or_voltage_limit)

            busNumberToDataOfBus[bus_number] = current_bus_data

            if (current_bus_data.Type == TypeOfBara.Slack):
                SlackBara = current_bus_data
            elif (current_bus_data.Type == TypeOfBara.PQ):
                PQBaras.append(current_bus_data)
            elif (current_bus_data.Type == TypeOfBara.PV):
                PVBaras.append(current_bus_data)

            # YBUS constructing for busses (constructing diagonals)
            i = busNumberToDataOfBus[bus_number].BusIndex

            YBUS[i, i] += np.complex128(bus_conductance + bus_susceptance * 1j)

            busCounter = busCounter + 1
            pass
        ### reading branch data ###
        elif (readingMode == StateOfReader.BranchReadingMode):
            # debugging the data
            # print("----------------------------------------")
            # print("-- i bus --")
            # print(line[0:4].strip())
            ibus = int(line[0:4].strip())
            # print(ibus)
            # print("-- j bus --")
            # print(line[5:9].strip())
            jbus = int(line[5:9].strip())
            # print(jbus)
            # print("--branch line charging--")
            # print(line[40:50].strip())
            branch_line_charging = float(line[40:50].strip())
            # print(branch_line_charging)
            # print("--branch reactance--")
            # print(line[29:40].strip())
            branch_reactance = float(line[29:40].strip())
            # print(branch_reactance)
            # print("--branch resistance--")
            # print(line[19:29].strip())
            branch_resistance = float(line[19:29].strip())
            # print(branch_resistance)
            # print("--transformer final turn ratio--")
            # print(line[76:82].strip())
            a_magnitude = float(line[76:82].strip())
            # print(a_magnitude)
            # print("--transformer (phase shifter)--")
            # print(line[83:90].strip())
            a_phase = float(line[83:90].strip()) * mat.pi / 180  # in radians
            # print(a_phase)
            # print("----------------------------------------")

            # if a_magnitude is 0 it is taken as 1
            if (a_magnitude < 1e-10):
                a_magnitude = 1

            # YBUS constructing for bus pairs (branches)
            i = busNumberToDataOfBus[ibus].BusIndex
            j = busNumberToDataOfBus[jbus].BusIndex

            a = a_magnitude * np.complex128(mat.cos(a_phase) + mat.sin(a_phase) * 1j)
            a_conj = np.conj(a)

            impedance = np.complex128(branch_resistance + branch_reactance * 1j)
            halfOfBranchLineSusceptance = np.complex128(branch_line_charging * 1j) / 2

            YBUS[i, i] += 1 / (impedance * a_magnitude * a_magnitude) + halfOfBranchLineSusceptance
            YBUS[i, j] += -1 / (impedance * a_conj)
            YBUS[j, i] += -1 / (impedance * a)
            YBUS[j, j] += 1 / impedance + halfOfBranchLineSusceptance

        previousLineDataHolderForGettingPowerBase = line

    magnitudes = np.abs(YBUS)

    marker_size = 1
    if(busCounter < 51):
        marker_size = 6
    elif(busCounter < 101):
        marker_size = 5
    elif(busCounter < 151):
        marker_size = 4
    elif(busCounter < 201):
        marker_size = 3
    elif(busCounter < 251):
        marker_size = 2


    # Plot the sparsity pattern
    plt.figure(figsize=(6,8))
    plt.spy(magnitudes, markersize = marker_size)
    plt.title('Sparsity Pattern of YBUS Magnitudes for ' + cdfDirectory[0:-4])
    plt.xlabel('Column Index')
    plt.ylabel('Row Index')
    plt.show()




    return YBUS, SlackBara, PVBaras, PQBaras








def Run(directory,MAXITERATION, EPSILON):
    # reading data and creating YBUS also the type of the baras
    (YBUS, SlackBara, PVBaras, PQBaras) = Read(directory)
    # running the Newton Raphsen method
    start_time = time.time()
    (currentVoltages, currentAngles, numberOfIterations) = Simulate(YBUS, SlackBara, PVBaras, PQBaras, MAXITERATION, EPSILON, False)
    end_time = time.time()
    solution_time = end_time - start_time
    # controlling the datas obtained from the simulation
    Control(currentVoltages, currentAngles, SlackBara, PVBaras, PQBaras)
    # calculating loss
    (realLoss, reactiveLoss) = CalculateLoss(SlackBara,PVBaras,PQBaras,currentVoltages,currentAngles,YBUS)
    # printing required data to console
    Print(currentVoltages, currentAngles, solution_time, numberOfIterations, realLoss, reactiveLoss, PVBaras, PQBaras, SlackBara)



if __name__ == "__main__":


    Run("ieee300cdf.txt", 100, 5e-13)



    # execution_time = end_time - start_time
    # print("Execution time:", execution_time, "seconds")


















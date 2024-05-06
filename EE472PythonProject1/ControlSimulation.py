
import numpy as np
import cmath as mat

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
    print("Indices where the difference exceeds 1e-2:")
    for index in zip(indices[0]):
        IsThereError = True
        print("index: " + str(index) + " calculated voltage: " + str(currentVoltages[index]) + " real voltage: " + str(
            TrueVoltages[index]))

    # Compute the absolute difference between the column vectors
    differenceAngles = np.abs(currentAngles - TrueAngles) * 180 / mat.pi

    # Find indices where the difference exceeds 1e-10
    indices = np.where(differenceAngles > 1e-0)

    # Print the indices where the difference exceeds the threshold voltages
    print("Indices where the difference exceeds 1e-0:")
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
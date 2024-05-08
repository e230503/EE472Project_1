import Constants as const
import numpy as np
import cmath as mat
import matplotlib.pyplot as plt

def Read(cdfDirectory):

    ieeecdf = open(cdfDirectory,"r")

    ########## creating YBus matrix ##########
    readingMode = const.StateOfReader.EmpthyReadingMode
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
            readingMode = const.StateOfReader.BranchReadingMode
            continue
        elif (line.__contains__("BUS DATA FOLLOWS ")):
            readingMode = const.StateOfReader.BusReadingMode

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
            readingMode = const.StateOfReader.EmpthyReadingMode
            continue

        ##### reading section #####
        ### reading bus data ###
        if (readingMode == const.StateOfReader.BusReadingMode):
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
            type_of_bara = const.TypeOfBara.Slack
            if (type_of_bara_int == 0 or type_of_bara_int == 1):
                type_of_bara = const.TypeOfBara.PQ
            elif (type_of_bara_int == 2):
                type_of_bara = const.TypeOfBara.PV
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

            current_bus_data = const.DataOfBusForSimulation(bus_number, busCounter, type_of_bara, final_voltage,
                                                            final_angle, load_MW, load_MVAR, generation_MW,
                                                            generation_MVAR, max_MVAR_or_voltage_limit,
                                                            min_MVAR_or_voltage_limit)

            busNumberToDataOfBus[bus_number] = current_bus_data

            if (current_bus_data.Type == const.TypeOfBara.Slack):
                SlackBara = current_bus_data
            elif (current_bus_data.Type == const.TypeOfBara.PQ):
                PQBaras.append(current_bus_data)
            elif (current_bus_data.Type == const.TypeOfBara.PV):
                PVBaras.append(current_bus_data)

            # YBUS constructing for busses (constructing diagonals)
            i = busNumberToDataOfBus[bus_number].BusIndex

            YBUS[i, i] += np.complex128(bus_conductance + bus_susceptance * 1j)

            busCounter = busCounter + 1
            pass
        ### reading branch data ###
        elif (readingMode == const.StateOfReader.BranchReadingMode):
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

    TotalP = 0
    TotalQ = 0

    TotalP = TotalP + SlackBara.GenerationP
    TotalQ = TotalQ + SlackBara.GenerationQ - SlackBara.LoadQ

    for PVbara in PVBaras:
        TotalP = TotalP + PVbara.GenerationP
        TotalQ = TotalQ + PVbara.GenerationQ - PVbara.LoadQ
    for PQBara in PQBaras:
        TotalP = TotalP - PQBara.LoadP
        TotalQ = TotalQ + PQBara.GenerationQ - PQBara.LoadQ

    print("-----------------------------------------------------")
    print("TotalP: " + str(TotalP) + " TotalQ: " + str(TotalQ))
    print("-----------------------------------------------------")


    return YBUS, SlackBara, PVBaras, PQBaras

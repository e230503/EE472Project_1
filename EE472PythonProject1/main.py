import DataReaderandYBUSCreater as R
import NRSimulation as S
import ControlSimulation as C
import time


def Run(directory,MAXITERATION, EPSILON):
    # reading data and creating YBUS also the type of the baras
    (YBUS, SlackBara, PVBaras, PQBaras) = R.Read(directory)
    # running the Newton Raphsen method
    (currentVoltages, currentAngles) = S.Simulate(YBUS, SlackBara, PVBaras, PQBaras, MAXITERATION, EPSILON, False)
    # controlling the datas obtained from the simulation
    C.Control(currentVoltages, currentAngles, SlackBara, PVBaras, PQBaras)


if __name__ == "__main__":
    start_time = time.time()

    Run("ieee57cdf.txt", 25, 5e-13)

    end_time = time.time()

    execution_time = end_time - start_time

    print("Execution time:", execution_time, "seconds")


















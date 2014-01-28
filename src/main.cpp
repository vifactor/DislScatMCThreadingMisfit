#include <mpi.h>
#include "RootEngine.h"

enum {rootTASK = 0};

static void master(std::string cfgname);
static void slave(std::string cfgname, int id);

int main(int argc, char** argv)
{
    int taskid;

    /*-------------mpi section-------------*/
	MPI_Init(&argc, &argv);

	MPI_Comm_rank(MPI_COMM_WORLD,&taskid);
    if (taskid == rootTASK)
    {
        master("default.cfg");
    }
    else
    {
        slave("default.cfg", taskid);
    }

    MPI_Finalize();
	/*-------------end of mpi section-------------*/
    return 0;
}

static void master(std::string filename)
{
    int ntasks;
    RootEngine engine;
    MCCalculator * calculator;
    MCCalculator::Data * curData, * cumData;
    int buffsize;
    double * recvbuff, *sendbuff;
    const double multiplicator = 1.2;

    /*allocations and initializations*/
    engine.setup(filename, 0);
    calculator = engine.allocateMCCalculator();
    curData = engine.allocateMCCalculatorData();
    cumData = engine.allocateMCCalculatorData();
    buffsize = 2 * curData->size();
    recvbuff = new double[buffsize];
    sendbuff = new double[buffsize];
    cumData->init();

    /* Find out how many processes there are in the default
     communicator */
    MPI_Comm_size(MPI_COMM_WORLD, &ntasks);

    /*wait while all preparations are finished*/
    MPI_Barrier(MPI_COMM_WORLD);

    /* Loop over getting new work requests until there is no more work
     to be done */
    while (!engine.done(cumData))
    {
        /* Send the slave a new work unit */
        MPI_Bcast(&curData->m_nbSteps, 1, MPI_INT, 0, MPI_COMM_WORLD);

        /*master also works*/
        calculator->run(curData);
        curData->transferTo(sendbuff);

        /* Receive results from slaves */
        MPI_Reduce(sendbuff, recvbuff, buffsize, MPI_DOUBLE, MPI_SUM,
					rootTASK, MPI_COMM_WORLD );

        /**append cummulative data and save result**/
        cumData->append(recvbuff, curData->m_nbSteps * ntasks);
        engine.save(cumData);

        /* increase the number of steps to perform */
        curData->m_nbSteps *= multiplicator;
    }
    /* Tell all the slaves to exit by sending 0 nb of steps to perform*/
    curData->m_nbSteps = 0;
    /* Send the slave a new work unit */
    MPI_Bcast(&curData->m_nbSteps, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);

    engine.freeMCCalculator(calculator);
    engine.freeMCCalculatorData(curData);
    delete[] recvbuff;
    delete[] sendbuff;

    std::cout << "Done." << std::endl;
}

static void slave(std::string filename, int taskid)
{
    Engine engine;
    MCCalculator * calculator;
    MCCalculator::Data * curData;
    int buffsize;
    double * recvbuff, *sendbuff;

    /*allocations and initializations*/
    engine.setup(filename, taskid);
    calculator = engine.allocateMCCalculator();
    curData = engine.allocateMCCalculatorData();
    buffsize = 2 * curData->size();
    recvbuff = new double[buffsize];
    sendbuff = new double[buffsize];

    MPI_Barrier(MPI_COMM_WORLD);
    while (true)
    {
        /* Receive a message from the master */
        MPI_Bcast(&curData->m_nbSteps, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);

        /* Check the tag of the received message. */
        if (curData->m_nbSteps == 0)
            break;

        /* Do the work */
        calculator->run(curData);
        curData->transferTo(sendbuff);

        /* Send the result back */
        MPI_Reduce(sendbuff, recvbuff, buffsize, MPI_DOUBLE, MPI_SUM,
					rootTASK, MPI_COMM_WORLD );
    }
    engine.freeMCCalculator(calculator);
    engine.freeMCCalculatorData(curData);
    delete[] recvbuff;
    delete[] sendbuff;
}

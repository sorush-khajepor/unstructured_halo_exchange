template <class T>
class Communication
{
private:
    Block2d<T> &block;
    const std::vector<Neighbor> neighbours;
    std::vector<MPI_Request> sendRequests;
    std::vector<MPI_Request> recvRequests;

public:
    Communication(Block2d<double> &block_, const std::vector<Neighbor> &neighbours_) : block{block_}, neighbours{neighbours_} {};

    auto Communicate(const Neighbor &nei)
    {
        MPI_Request recvRequest;
        MPI_Request sendRequest;

        MPI_Irecv(&block(0, 0), 1, nei.inType, nei.rank, nei.inTag, MPI_COMM_WORLD, &recvRequest);
        MPI_Isend(&block(0, 0), 1, nei.outType, nei.rank, nei.outTag, MPI_COMM_WORLD, &sendRequest);
        return std::array<MPI_Request, 2>{sendRequest, recvRequest};
    }

    void CommunicateAsync()
    {
        for (auto &neighbour : neighbours)
        {
            auto [sendReq, recvReq] = Communicate(neighbour);
            sendRequests.emplace_back(sendReq);
            recvRequests.emplace_back(recvReq);
        }
    }

    void AwaitRecv()
    {
        for (auto &req : recvRequests)
            MPI_Wait(&req, MPI_STATUS_IGNORE);
        recvRequests.clear();
    }

    void AwaitSend()
    {
        for (auto &req : sendRequests)
            MPI_Wait(&req, MPI_STATUS_IGNORE);
        sendRequests.clear();
    }
    void Await(){
        AwaitRecv();
        AwaitSend();
    }

    void Communicate()
    {
        CommunicateAsync();
        Await();
        
    }
};

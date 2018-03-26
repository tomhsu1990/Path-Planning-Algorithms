# Path-Planning-Algorithms

Path planning (also known as the navigation problem or the piano mover's problem) is a fundemental problem in Robotics to break down a desired movement task into discrete motions that satisfy movement constraints and possibly optimize some aspect of the movement. There are three main approaches in Path Planning:

1. The Exact approach which can, in theory, solve any algebraic planning problem. However, practitioners tend to implement exact algorithms using machine approximations; then it is no longer clear what the guarantees of exact algorithms mean. In short, "Exact" algorithms are impractical except for the simplest cases.

2. The Sampling approach has been dominant in path planning research domain for the last two decades. We refer to Probabilistic Road Map (PRM) and Rapidly-exploring Random Tree (RRT) the two main representatives. The greatest strength is that it can deal with configuration in high dimension and is relatively easy to implement. But, Its central problem is the inability to halt (couched as “narrow passage problem”).

3. The Cell Decomposition approach is the oldest among the three approaches, and remains popular with practitioners. The basic idea behind this method is that a path between the initial configuration and the goal configuration can be determined by subdividing the free space of the robot's configuration into smaller regions called cells. From the constructed connectivity graph (cells), a continuous path, or channel, can be determined by simply following adjacent free cells from the initial point to the goal point. The main advantage is that it guarantees the completeness. That is, if there exists a free path, the Cell Decomposition approach will find it.

from math import inf
import copy
import sys
from math import radians, cos, sin, sqrt, atan2

# Student Name: Samuel Emmett

n = 0
class WeightedAdjacencyMatrix :
    """A weighted graph represented as a matrix."""

    __slots__ = ['_W']

    def __init__(self, size, edges=[], weights=[]) :
        """Initializes a weighted adjacency matrix for a graph with size nodes.

        Graph is initialized with size nodes and a specified set of
        edges and edge weights.
        
        Keyword arguments:
        size -- Number of nodes of the graph.
        edges -- a list of ordered pairs (2-tuples) indicating the
                 edges of the graph.  The default value is an empty list
                 which means no edges by default.
        weights -- a list of weights for the edges, which should be the same
                   length as the edges list.  The position of a value in
                   the weights list corresponds to the edge in the same
                   position of the edges list.
        """
        
        self._W = [ [] for i in range (size) ] #initializes a completely empty list of n size

        #This initialized the empty weighted matrix
        for x in range(size):
            for y in range(size):
                if x == y: #if x and y are = then the distance from x to y is 0
                    self._W[x].append(0)
                else : #Otherwise the distance from x to y will be initialized as infinity
                    self._W[x].append(inf)

        #Adds edges and weights to the matrix
        for i in range (len(weights)): #since weights and edges need to be of equal length it doesn't matter which is used here
           u,v = edges[i]
           self.addEdge(u,v,weights[i])


    def addEdge(self, u, v, weight) :
        """Adds an undirected edge between u to v with the specified weight.

        Keyword arguments:
        u -- vertex id (0-based index)
        v -- vertex id (0-based index)
        weight -- edge weight
        """

        #You will need to add matrix[u][v] and matrix[v][u]

        self._W[u][v] = weight
        self._W[v][u] = weight

    def floydWarshall(self) :
        """Floyd Warshall algorithm for all pairs shortest paths.

        Returns a matrix D consisting of the weights of the shortest
        paths between all pairs of vertices, and a matrix P for
        the predecessors matrix (what the textbook called PI).
        This method MUST NOT change the weight matrix of the graph
        itself.  
        """
        # Your return statement will look something like this one
        # in the comment on the following line.  That is return
        # the two matrices, with the D matrix first.  The return None
        # is just a placeholder so that this is valid Python syntax before
        # you've completed the assignment.  This comment line is
        # more like what it should look like:
        # return D, P

        n = len(self._W)
        D = copy.deepcopy(self._W) #this way self._W will not be altered by the method
        P = [ [ None for i in range(n) ] for j in range(n) ]

        for i in range (n) :
            for j in range (n) :
                if i != j and self._W[i][j] < inf :
                    P[i][j] = i
                else :
                    P[i][j] = None

        for k in range (n) :
            for i in range (n) :
                for j in range (n) :
                    if D[i][k] + D[k][j] < D[i][j] :
                        D[i][j] = D[i][k] + D[k][j]
                        P[i][j] = P[k][j]

        return D, P

class WeightedDirectedAdjacencyMatrix(WeightedAdjacencyMatrix) :
    """A weighted digraph represented as a matrix."""

    def addEdge(self, u, v, weight) :
        """Adds a directed edge from u to v with the specified weight.

        Keyword arguments:
        u -- source vertex id (0-based index)
        v -- target vertex id (0-based index)
        weight -- edge weight
        """
        self._W[u][v] = weight

def testFloydWarshall() :
    """See assignment instructions at top."""
    weights = [2, 4, 1, 2, 3, 2]
    edges = [(1, 0), (0, 2), (1, 3), (1, 4), (2, 3), (3, 4)]

    graph = WeightedAdjacencyMatrix(5, edges, weights)
    D1, P1 = graph.floydWarshall()

    print("FOR AN UNDIRECTED MATRIX")

    print("The Pi Matrix is:")
    print('\n'.join([''.join(['{:^5}'.format(str(item)) for item in row]) for row in P1]))
    print()
    print("The D Matrix is:")
    print('\n'.join([''.join(['{:^5}'.format(str(item)) for item in row]) for row in D1]))
    print()

    graph2 = WeightedDirectedAdjacencyMatrix(5, edges, weights)
    D2, P2 = graph2.floydWarshall()

    print("FOR A DIRECTED MATRIX")

    print("The Pi Matrix is:")
    print('\n'.join([''.join(['{:^5}'.format(str(item)) for item in row]) for row in P2]))
    print()
    print("The D Matrix is:")
    print('\n'.join([''.join(['{:^5}'.format(str(item)) for item in row]) for row in D2]))
    print()

def parseHighwayGraphMatrix(filename) :
    """EXTRA CREDIT: Rewrite your highway graph parser from
    assignment 2 here in this function but return a WeightedAdjacencyMatrix
    object from this function. If you had that assignment working,
    then this part of the extra credit should be very easy (i.e.,
    copying and pasting code and then making very minor adjustments
    to use construct and return a WeightedAdjacencyMatrix object
    instead of the other graph type you already have.

    Keyword arguments:
    filename - the name of a highhway graph file.
    """

    def haversine(lat1, lng1, lat2, lng2):
        """Computes haversine distance between two points in latitude, longitude.

        Keyword Arguments:
        lat1 -- latitude of point 1
        lng1 -- longitude of point 1
        lat2 -- latitude of point 2
        lng2 -- longitude of point 2

        Returns haversine distance in meters.
        """
        R = 6371e3  # radius of the earth in meters
        phi1 = radians(lat1)
        phi2 = radians(lat2)
        deltaPhi = radians((lat2 - lat1))
        deltaLam = radians((lng2 - lng1))

        a = (sin(deltaPhi / 2) * sin(deltaPhi / 2)) + ((cos(phi1) * cos(phi2)) * sin(deltaLam / 2) * sin(deltaLam / 2))
        c = 2 * atan2(sqrt(a), sqrt(1 - a))
        d = R * c
        return d  # Returns the distance in meters

    file = open(filename)
    count = 0
    vertex = []

    weights = []
    edgeData = []

    while True:  # READS LINES
        count += 1
        line = file.readline()
        if line == "":
            file.close()
            break
        if count == 2:  # uses second line for vertex and edge numbers
            vertexNum, edgeNum = line.split()
            size = int(vertexNum)
            edgeSize = int(edgeNum)

            # May need to be offset
        if count >= 3:  # needed or else the statement will run on the first line and second line

            # This reads the latitude and longitude data
            if (count <= size + 2) and (count != size + 3):
                a0, a1, a2 = line.split()
                v = (float(a1), float(a2))
                vertex.append(v)

                # Maybe Add the verticese data to the verticese list, then read the data at the list[p] section
                # where p = the specified vertices by a(o) and (a1)

            if count > size + 2 and line != "":
                a0, a1, a2 = line.split()
                edge = (int(a0), int(a1))
                lat1, lng1 = vertex[int(a0)]
                lat2, lng2 = vertex[int(a1)]
                weights.append(haversine(lat1, lng1, lat2, lng2))
                edgeData.append(edge)

    HighGraph = WeightedAdjacencyMatrix(size, edgeData, weights)

    return HighGraph

def pairShortestPath(D, P, s, t) :
    """EXTRA CREDIT: This function takes D and P matrices (i.e., what is generated
    by floydWarshall), and a source vertex (where you want to start) and
    destination or target vertex t (where you want to end up) and
    returns a pair: w, path, such that w is the weight of the shortest
    path from s to t (just a simple lookup in the D matrix) and
    path is a Python list of vertex ids starting at s and ending at t
    derived from the P matrix.

    Keyword arguments:
    D - the D matrix
    P - the Pi matrix
    s - the source vertex
    t - the destination vertex
    """
    # Your actual return will look something like this:
    # return w, path

    row = int(s)
    col = int(t)
    w = D[row][col]
    path = []

    if w == None :
        return w, path
    for x in range (len(P)): #at most the loop will run O(x) times
        if P[row][col] == int(s) :
            path.insert(0, int(s))
            path.append(int(t))
            break
        elif P[row][col] != None :
            path.insert(0, P[row][col])
            col = P[row][col]

    return w, path

if __name__ == "__main__" :
    # EXTRA CREDIT: Write code here that:
    #   (a) Gets the name of a highway graph file from the command line
    #       arguments.
    #   (b) Uses parseHighwayGraphMatrix from above to parse that file
    #       into a WeightedAdjacencyMatrix object.
    #   (c) Runs the floydWarshall method on that graph.
    #   (d) Then, prompts the user (use the Python docs to figure out
    #       how to do this) for a source and target vertex, s and t.
    #   (e) Uses pairShortestPath to get the weight of the shortest
    #       path between their chosen s and t, and the path itself.
    #   (f) Outputs the weight and path.
    #   (g) Repeats d, e, and f in a loop until the user indicates they
    #       want to quit.  You can decide how to get that decision from them.
                #0      1      2    3   4
    i = 0
    for arg in sys.argv:
        if i != 0:
            HighwayGraph = parseHighwayGraphMatrix(str(sys.argv[i]))
            D, P = HighwayGraph.floydWarshall()
            print("The Pi Matrix of the HighwayMatrix is:")
            print('\n'.join([''.join(['{:^11}'.format(str(item)[0:10]) for item in row]) for row in P]))
            print("")
            print("The D Matrix of the HighwayMatrix is:")
            print('\n'.join([''.join(['{:^11}'.format(str(item)[0:10]) for item in row]) for row in D]))
            print("Please be aware that this Matrix will trunkate some values ")
            print("")
            print("")
            ans = "y"

            while True:
                print("What vertex would you like to start at? Please choose a number between 0 and " + str(len(D)-1))
                s = input("Where would you like to start?:")
                while not s.isnumeric() :
                    s = input("Please input an integer in the specified range:")
                print()

                print("What vertex would you like to traverse to? Please choose a number between 0 and " + str(len(D)-1))
                t = input("Where would you like to go?:")
                while not t.isnumeric() : #Makes sure that the input is numeric
                    t = input("Please input an integer in the specified range:")
                print()

                print(pairShortestPath(D, P, s, t))

                ans = input("Would you like to run the program again? (y/n) [case-sensitive]:")
                while ans != "y" and ans != "n": #Makes sure that the input is numeric
                    ans = input("Please type y or n [case-sensitive]:")
                if ans != "y":
                    break
            i = i + 1
        else:
            i = i + 1



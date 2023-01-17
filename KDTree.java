
/**
* JohnC 2023
*/
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.function.BiFunction;
import java.util.function.Function;

/**
 * this is the KDTree class
 * 
 * @author John C
 */
class KDTree {
    KDTree self;
    double[][] points;
    BiFunction<double[], double[], Double> metric;
    int dimensions;
    Node root;

    /**
     * init a instance with just points
     * 
     * @param points a 2D array like this[p1, p2, ...] with p1: [x1, y1, z1, ...]
     */
    public KDTree(double[][] points) {
        this(points, null, 2);
    }

    /**
     * init a instance with more options
     * 
     * @param points     a 2D array like this[p1, p2, ...] with p1: [x1, y1, z1,
     *                   ...]
     * @param metric     the metric function you wanna use for the tree
     * @param dimensions dimenstions of your data, 2, 3, etc.
     */
    public KDTree(double[][] points, BiFunction<double[], double[], Double> metric, int dimensions) {
        this.self = this;
        this.points = points;
        this.metric = metric;
        if (this.metric == null) {
            this.metric = (pointA, pointB) -> (pointA[0] - pointB[0]) * (pointA[0] - pointB[0])
                    + (pointA[1] - pointB[1]) * (pointA[1] - pointB[1]);
        }
        this.dimensions = dimensions;
        this.root = this.buildTree(points, 0, null);
    }

    /**
     * build the kd tree
     * 
     * @param points data points, a 2D array like this[p1, p2, ...] with p1: [x1,
     *               y1, z1, ...]
     * @param depth  the depth count of the tree
     * @param parent the root of the new tree
     * @return the root node of the tree
     */
    public Node buildTree(double[][] points, int depth, Node parent) {
        int dim = depth % dimensions;
        int median;
        Node node;
        if (points == null || points.length == 0) {
            return null;
        }
        if (points.length == 1) {
            return new Node(points[0], dim, parent);
        }

        Arrays.sort(points, (a, b) -> {
            return (int) (a[dim] - b[dim]);
        });

        median = (int) Math.floor(points.length / 2);
        node = new Node(points[median], dim, parent);
        node.left = this.buildTree(Arrays.copyOfRange(points, 0, median), depth + 1, node);
        node.right = this.buildTree(Arrays.copyOfRange(points, median + 1, points.length), depth + 1, node);

        return node;
    }

    private Node _innerSearch(double[] point, Node node, Node parent) {
        if (node == null) {
            return parent;
        }
        int dimension = node.dimension;
        if (point[dimension] < node.obj[dimension]) {
            return this._innerSearch(point, node.left, node);
        } else {
            return this._innerSearch(point, node.right, node);
        }
    }

    /**
     * insert a point to the tree
     * 
     * @param point point [x, y, z, ...] to be inserted
     */
    public void insert(double[] point) {
        Node insertPosition = this._innerSearch(point, this.root, null);
        Node newNode;
        int dimension;

        if (insertPosition == null) {
            this.root = new Node(point, 0, null);
            return;
        }

        newNode = new Node(point, (insertPosition.dimension + 1) % this.dimensions, insertPosition);
        dimension = insertPosition.dimension;

        if (point[dimension] < insertPosition.obj[dimension]) {
            insertPosition.left = newNode;
        } else {
            insertPosition.right = newNode;
        }
    }

    private Node _nodeSearch(double[] point, Node node) {
        if (node == null)
            return null;
        if (Arrays.equals(node.obj, point)) {
            return node;
        }
        int dimension = node.dimension;
        if (point[dimension] < node.obj[dimension]) {
            return this._nodeSearch(point, node.left);
        } else {
            return this._nodeSearch(point, node.right);
        }
    }

    private Node _findMin(Node node, int dim) {
        int dimension;
        double own;
        Node left;
        Node right;
        Node min;
        if (node == null)
            return null;
        dimension = dim;
        if (node.dimension == dim) {
            if (node.left != null) {
                return this._findMin(node.left, dim);
            }
            return node;
        }

        own = node.obj[dimension];
        left = this._findMin(node.left, dim);
        right = this._findMin(node.right, dim);
        min = node;
        if (left != null && left.obj[dimension] < own) {
            min = left;
        }
        if (right != null && right.obj[dimension] < min.obj[dimension]) {
            min = right;
        }
        return min;
    }

    private void _removeNode(double[] point, Node node) {
        Node nextNode;
        double[] nextObj;
        int pDimension;
        if (node.left == null && node.right == null) {
            if (node.parent == null) {
                this.self.root = null;
                return;
            }

            pDimension = node.parent.dimension;

            if (node.obj[pDimension] < node.parent.obj[pDimension]) {
                node.parent.left = null;
            } else {
                node.parent.right = null;
            }
            return;
        }
        if (node.right != null) {
            nextNode = this._findMin(node.right, node.dimension);
            nextObj = nextNode.obj;
            this._removeNode(point, nextNode);
            node.obj = nextObj;
        } else {
            nextNode = this._findMin(node.left, node.dimension);
            nextObj = nextNode.obj;
            this._removeNode(point, nextNode);
            node.right = node.left;
            node.left = null;
            node.obj = nextObj;
        }

    }

    /**
     * remove a point from the tree
     * 
     * @param point point [x, y, z, ...] to remove
     */
    public void remove(double[] point) {
        Node node;
        node = this._nodeSearch(point, this.root);

        if (node == null)
            return;

        this._removeNode(point, node);
    }

    private void _saveNode(BinaryHeap bestNodes, int maxNodes, Node node, double distance) {
        HashMap<String, Object> tem = new HashMap<String, Object>();
        tem.put("0", node);
        tem.put("1", distance);
        bestNodes.push(tem);
        if (bestNodes.size() > maxNodes) {
            bestNodes.pop();
        }
    }

    private void _nearestSearch(double[] point, int maxNodes, BinaryHeap bestNodes, Node node,
            BiFunction<double[], double[], Double> metric, int dimensions) {
        Node bestChild;
        int dimension = node.dimension;
        double ownDistance = metric.apply(point, node.obj).doubleValue();
        double[] linearPoint = new double[dimensions];
        double linearDistance;
        Node otherChild;
        int i;

        for (i = 0; i < dimensions; i++) {
            if (i == node.dimension) {
                linearPoint[i] = point[i];
            } else {
                linearPoint[i] = node.obj[i];
            }
        }

        linearDistance = metric.apply(linearPoint, node.obj);
        if (node.right == null && node.left == null) {
            if (bestNodes.size() < maxNodes || ownDistance < (Double) bestNodes.peek().get("1")) {
                this._saveNode(bestNodes, maxNodes, node, ownDistance);
            }
            return;
        }

        if (node.right == null) {
            bestChild = node.left;
        } else if (node.left == null) {
            bestChild = node.right;
        } else {
            if (point[dimension] < node.obj[dimension]) {
                bestChild = node.left;
            } else {
                bestChild = node.right;
            }
        }

        this._nearestSearch(point, maxNodes, bestNodes, bestChild, metric, dimensions);

        if (bestNodes.size() < maxNodes || ownDistance < (Double) bestNodes.peek().get("1")) {
            this._saveNode(bestNodes, maxNodes, node, ownDistance);
        }

        if (bestNodes.size() < maxNodes || Math.abs(linearDistance) < (Double) bestNodes.peek().get("1")) {
            if (bestChild == node.left) {
                otherChild = node.right;
            } else {
                otherChild = node.left;
            }
            if (otherChild != null) {
                this._nearestSearch(point, maxNodes, bestNodes, otherChild, metric, dimensions);
            }
        }
    }

    @SuppressWarnings("unchecked")
    /**
     * return nearest point(s) of a point
     * 
     * @param point       the point
     * @param maxNodes    max amount of neibours returned
     * @param maxDistance max distance of the neibours from the point returned
     * @return an arry of HashMap, Each map: "obj"(double[]) : the neibour point;
     *         "dist"(double) : the distance from this point to the input point
     */
    public HashMap<String, Object>[] nearest(double[] point, int maxNodes, double maxDistance) {
        int i;
        ArrayList<HashMap<String, Object>> result = new ArrayList<HashMap<String, Object>>();
        BinaryHeap bestNodes = new BinaryHeap((e) -> -((Double) e.get("1")).doubleValue());

        if (maxDistance > 0) {
            for (i = 0; i < maxNodes; i++) {
                HashMap<String, Object> tem = new HashMap<String, Object>();
                tem.put("0", null);
                tem.put("1", maxNodes);
                bestNodes.push(tem);
            }
        }

        if (this.root != null)
            this._nearestSearch(point, maxNodes, bestNodes, this.root, this.metric, this.dimensions);

        for (i = 0; i < Math.min(maxNodes, bestNodes.content.size()); i++) {
            if (bestNodes.content.get(i).get("0") != null) {
                Node t = (Node) bestNodes.content.get(i).get("0");
                HashMap<String, Object> tem = new HashMap<String, Object>();
                tem.put("obj", t.obj);
                tem.put("dist", (Double) bestNodes.content.get(i).get("1"));
                result.add(tem);
            }
        }

        return result.toArray((HashMap<String, Object>[]) new HashMap[result.size()]);
    }
}

class Node {
    public double[] obj;
    public Node left;
    public Node right;
    public Node parent;
    int dimension;

    public Node(double[] obj, int dimension, Node parent) {
        this.obj = obj;
        this.left = null;
        this.right = null;
        this.parent = parent;
        this.dimension = dimension;
    }
}

class BinaryHeap {
    public ArrayList<HashMap<String, Object>> content;
    public Function<HashMap<String, Object>, Double> scoreFunction;

    public BinaryHeap(Function<HashMap<String, Object>, Double> scoreFunction) {
        this.content = new ArrayList<HashMap<String, Object>>();
        this.scoreFunction = scoreFunction;
    }

    public void push(HashMap<String, Object> element) {
        this.content.add(element);
        this.bubbleUp(this.content.size() - 1);
    }

    public HashMap<String, Object> pop() {
        HashMap<String, Object> result = this.content.get(0);
        HashMap<String, Object> end = this.content.remove(this.content.size() - 1);
        if (this.content.size() > 0) {
            this.content.set(0, end);
            this.sinkDown(0);
        }
        return result;
    }

    public HashMap<String, Object> peek() {
        return this.content.get(0);
    }

    public void remove(int idx) {
        int len = this.content.size();
        for (int i = 0; i < len; i++) {
            if (i == idx) {
                HashMap<String, Object> end = this.content.remove(this.content.size() - 1);
                HashMap<String, Object> node = this.content.get(i);
                if (i != len - 1) {
                    this.content.set(i, end);
                    if (this.scoreFunction.apply(end) < this.scoreFunction.apply(node)) {
                        this.bubbleUp(i);
                    } else {
                        this.sinkDown(i);
                    }
                }
                return;
            }
        }
        throw new RuntimeException("[Binary heap]: remove: node not found");
    }

    public int size() {
        return this.content.size();
    }

    public void bubbleUp(int n) {
        HashMap<String, Object> element = this.content.get(n);
        while (n > 0) {
            int parentN = (int) Math.floor((n + 1) / 2) - 1;
            HashMap<String, Object> parent = this.content.get(parentN);
            if (this.scoreFunction.apply(element) < this.scoreFunction.apply(parent)) {
                this.content.set(parentN, element);
                this.content.set(n, parent);
                n = parentN;
            } else {
                break;
            }
        }
    }

    public void sinkDown(int n) {
        int length = this.content.size();
        HashMap<String, Object> element = this.content.get(n);
        double elemScore = this.scoreFunction.apply(element);
        while (true) {
            int child2N = (n + 1) * 2;
            int child1N = child2N - 1;
            int swap = -100;
            double child1Score = Double.MIN_VALUE;
            if (child1N < length) {
                HashMap<String, Object> child1 = this.content.get(child1N);
                child1Score = this.scoreFunction.apply(child1);
                if (child1Score < elemScore)
                    swap = child1N;
            }
            if (child2N < length) {
                HashMap<String, Object> child2 = this.content.get(child2N);
                double child2Score = this.scoreFunction.apply(child2);
                if (child2Score < (swap < 0 ? elemScore : child1Score))
                    swap = child2N;
            }
            if (swap > 0) {
                this.content.set(n, this.content.get(swap));
                this.content.set(swap, element);
                n = swap;
            } else {
                break;
            }
        }
    }
}
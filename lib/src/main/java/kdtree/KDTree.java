package kdtree;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import java.util.function.BiFunction;
import java.util.function.Function;

/**
* JohnC 2024
for customized object with float as coordinate on each dimension, the class of that object need to have a nullary constructor
base on https://github.com/ubilabs/kd-tree-javascript/blob/master/kdTree.js
*/


public class KDTree {
    private float getDimension(Object point, String dimension) {
        try {
            return point.getClass().getField(dimension).getFloat(point);
        } catch (Exception e) {
            try {
                return point.getClass().getDeclaredField(dimension).getFloat(point);
            } catch (Exception ee) {
                throw new Error("bad dimension: " + dimension + "; object: " + point + "exception: " + ee);
            }
        }
    }

    private void setDimension(Object point, String dimension, float value) {
        try {
            point.getClass().getField(dimension).setFloat(point, value);
        } catch (Exception e) {
            try {
                point.getClass().getDeclaredField(dimension).setFloat(point, value);
            } catch (Exception ee) {
                throw new Error("bad dimension: " + dimension + "; object: " + point + "exception: " + ee);
            }
        }
    }

    private Object instantiate(Object point){
        if (processingSketch == null) {
            try {
                return point.getClass().getConstructor().newInstance();
            } catch (Exception e) {
                System.out.println(e);
                try {
                    return point.getClass().getDeclaredConstructor().newInstance();
                } catch (Exception ee) {
                   throw new Error("can't create new instance; class: " + point.getClass() + "; " + e);
                }
            }
        } else {
            try {
                return point.getClass().getConstructor(processingSketch.getClass()).newInstance(processingSketch);
            } catch (Exception e) {
                System.out.println(e);
                try {
                    return point.getClass().getDeclaredConstructor(processingSketch.getClass()).newInstance(processingSketch);
                } catch (Exception ee) {
                   throw new Error("can't create new instance; class: " + point.getClass() + "; " + e);
                }
            }
        }
    }

    private class Node {

        Object obj;
        int dimension;
        Node left, right;
        Node parent;

        Node(Object obj, int dimension, Node parent) {
            this.obj = obj;
            this.dimension = dimension;
            this.parent = parent;
            this.left = null;
            this.right = null;
        }
    }

    private class BinaryHeap {
        ArrayList<Map<String, Object>> content;
        Function<Map<String, Object>, Float> scoreFunction;

        BinaryHeap(Function<Map<String, Object>, Float> scoreFunction) {
            this.scoreFunction = scoreFunction;
            this.content = new ArrayList<>();
        }

        void push(Map<String, Object> element) {
            content.add(element);
            bubbleUp(content.size() - 1);
        }
        
        Map<String, Object> pop() {
            Map<String, Object> res = content.get(0);
            Map<String, Object> end = content.remove(content.size() - 1);
            if (content.size() > 0) {
                content.set(0, end);
                sinkDown(0);
            }
            return res;
        }

        int size() {
            return content.size();
        }

        Map<String, Object> peek() {
            return content.get(0);
        }

        void remove(Map<String, Object> node) {
            int len = content.size();
            for (int i = 0; i < len; i++) {
                if (content.get(i).equals(node)) {
                    Map<String, Object> end = content.remove(content.size() - 1);
                    if (i != len - 1) {
                        content.set(i, end);
                        if (scoreFunction.apply(end) < scoreFunction.apply(node)) {
                            bubbleUp(i);
                        } else {
                            sinkDown(i);
                        }
                    }
                    return;
                }
            }
            throw new Error("node Not found");
        }

        void bubbleUp(int n) {
            Map<String, Object> element = content.get(n);

            while (n > 0) {
                int parentN = (int) Math.floor((n + 1) / 2) - 1;
                Map<String, Object> parent = content.get(parentN);
                if (scoreFunction.apply(element) < scoreFunction.apply(parent)) {
                    content.set(parentN, element);
                    content.set(n, parent);
                    n = parentN;
                } else {
                    break;
                }
            }
        }

        void sinkDown(int n) {
            int length = content.size();
            Map<String, Object> element = content.get(n);
            float eleScore = scoreFunction.apply(element);
            float child1Score = 0, child2Score = 0;
            while (true) {
                int child2N = (n + 1) * 2, child1N = child2N - 1;
                int swap = -100;
                if (child1N < length) {
                    Map<String, Object> child1 = content.get(child1N);
                    child1Score = scoreFunction.apply(child1);
                    if (child1Score < eleScore) {
                        swap = child1N;
                    }
                }
                if (child2N < length) {
                    Map<String, Object> child2 = content.get(child2N);
                    child2Score = scoreFunction.apply(child2);
                    if (child2Score < (swap == -100 ? eleScore : child1Score)) {
                        swap = child2N;
                    }
                }
                if (swap != -100) {
                    content.set(n, content.get(swap));
                    content.set(swap, element);
                    n = swap;
                } else {
                    break;
                }
            }
        }
    }

    public String[] dimensions;
    public BiFunction<Object, Object, Float> metric;
    public Node root;
    public Object processingSketch;

    public KDTree(Object[] points, BiFunction<Object, Object, Float> metric, String[] dimensions, Object ps) {
        this.dimensions = dimensions;
        this.metric = metric;
        if (metric == null) {
            this.metric = (Object a, Object b) -> {
                float sum = 0;
                for (int i = 0; i < this.dimensions.length; i++) {
                    String fn = this.dimensions[i];
                    float vA = getDimension(a, fn);
                    float vB = getDimension(b, fn);
                    sum += (vA - vB) * (vA - vB);
                }
                return sum;
            };
        }
        this.root = buildTree(points, 0, null);
        this.processingSketch = ps;
    }

    public KDTree(Object[] points, BiFunction<Object, Object, Float> metric, String[] dimensions){
        this(points, metric, dimensions, null);
    }

    public Node buildTree(Object[] points, int depth, Node parent) {
        int dim = depth % this.dimensions.length;
        int median;
        Node node;
        if (points.length == 0)
            return null;
        if (points.length == 1)
            return new Node(points[0], dim, parent);
        Arrays.sort(points, (Object a, Object b) -> {
            String fn = this.dimensions[dim];
            float res;
            float vA = getDimension(a, fn);
            float vB = getDimension(b, fn);
            res = vA - vB;
            return res > 0 ? 1 : -1;
        });
        median = (int) Math.floor(points.length / 2);
        node = new Node(points[median], dim, parent);
        node.left = buildTree(Arrays.copyOfRange(points, 0, median), depth + 1, node);
        node.right = buildTree(Arrays.copyOfRange(points, median + 1, points.length), depth + 1, node);

        return node;
    }
    
    public void insert(Object point) {
        Node insertPosition = innerSearch(point, root, null);
        Node newNode;
        String dimension;
        if (insertPosition == null) {
            root = new Node(point, 0, null);
            return;
        }

        newNode = new Node(point, (insertPosition.dimension + 1) % dimensions.length, insertPosition);
        dimension = dimensions[insertPosition.dimension];
        float pV = getDimension(point, dimension), nV = getDimension(insertPosition.obj, dimension);
        if (pV < nV) {
            insertPosition.left = newNode;
        } else {
            insertPosition.right = newNode;
        }
    }
    
    private Node innerSearch(Object point, Node node, Node parent) {
        if (node == null)
            return parent;
        String dimension = this.dimensions[node.dimension];
        float pV = getDimension(point, dimension), nV = getDimension(node.obj, dimension);
        if (pV < nV) {
            return innerSearch(point, node.left, node);
        } else {
            return innerSearch(point, node.right, node);
        }
    }
    
    public void remove(Object point) {
        Node node;
        node = nodeSearch(point, root);
        if (node == null) return;
        removeNode(node);
    }

    private Node nodeSearch(Object point, Node node) {
        if (node == null) {
            return null;
        }
        if (node.obj.equals(point)) {
            return node;
        }

        String dimension = dimensions[node.dimension];
        float pV = getDimension(point, dimension), nV = getDimension(node.obj, dimension);
        if (pV < nV) {
            return nodeSearch(point, node.left);
        } else {
            return nodeSearch(point, node.right);
        }
    }

    private Node findMin(Node node, int dim) {
        String dimension;
        float own;
        Node left, right, min;
        if (node == null) return null;
        dimension = dimensions[dim];
        if (node.dimension == dim) {
            if (node.left != null)
                return findMin(node.left, dim);
            return node;
        }
        
        own = getDimension(node.obj, dimension);

        left = findMin(node.left, dim);
        right = findMin(node.right, dim);
        min = node;
        float minV = own;
        if (left != null) {
            float lV = getDimension(left.obj, dimension);
            if (lV < own) {
                min = left;
                minV = lV;
            }
        }
        if (right != null) {
            float rV = getDimension(right.obj, dimension);
            if (rV < minV) {
                min = right;
            }
        }
        return min;
    }

    private void removeNode(Node node) {
        Node nextNode;
        Object nextObject;
        String pDimension;
        if (node.left == null && node.right == null) {
            if (node.parent == null) {
                root = null;
                return;
            }

            pDimension = dimensions[node.parent.dimension];
            float nV = getDimension(node.obj, pDimension), pV = getDimension(node.parent.obj, pDimension);
            if (nV < pV) {
                node.parent.left = null;
            } else {
                node.parent.right = null;
            }
            return;
        }
        if (node.right != null) {
            nextNode = findMin(node.right, node.dimension);
            nextObject = nextNode.obj;
            removeNode(nextNode);
            node.obj = nextObject;
        } else {
            nextNode = findMin(node.left, node.dimension);
            nextObject = nextNode.obj;
            removeNode(nextNode);
            node.right = node.left;
            node.left = null;
            node.obj = nextObject;
        }
    }

    @SuppressWarnings("unchecked")
    public Map<String, Object>[] nearest(Object point, int maxNodes, float maxDist) {
        int i;
        ArrayList<Map<String, Object>> result = new ArrayList<>();
        BinaryHeap bestNodes;
        bestNodes = new BinaryHeap((e) -> -1 * (Float) e.get("dist"));
        if (maxDist > 0) {
            for (i = 0; i < maxNodes; i++) {
                Map<String, Object> tem = new HashMap<String, Object>();
                tem.put("obj", null);
                tem.put("dist", maxDist);
                bestNodes.push(tem);
            }
        }
        if (root != null) {
            nearestSearch(point, root, bestNodes, maxNodes);
        }

        for (i = 0; i < Math.min(maxNodes, bestNodes.content.size()); i++) {
            if (bestNodes.content.get(i).get("obj") != null) {
                Node n = (Node) bestNodes.content.get(i).get("obj");
                Map<String, Object> tem = new HashMap<String, Object>();
                tem.put("obj", n.obj);
                tem.put("dist", bestNodes.content.get(i).get("dist"));
                result.add(tem);
            }
        }
        
        return (Map<String, Object>[]) result.toArray(new Map[0]);
    }

    private void nearestSearch(Object point, Node node, BinaryHeap bestNodes, int maxNodes) {
        Node bestChild, otherChild;
        String dimension = dimensions[node.dimension];
        float ownDist = metric.apply(point, node.obj);
        Object linearPoint = instantiate(point);
        
        float linearDist;

        for (int i = 0; i < dimensions.length; i++) {
            if (i == node.dimension) {
                setDimension(linearPoint, dimensions[i], getDimension(point, dimensions[i]));
            } else {
                setDimension(linearPoint, dimensions[i], getDimension(node.obj, dimensions[i]));
            }
        }

        linearDist = metric.apply(linearPoint, node.obj);
        if (node.right == null && node.left == null) {
            if (bestNodes.size() < maxNodes || ownDist < (Float) bestNodes.peek().get("dist")) {
                saveNode(bestNodes, maxNodes, node, ownDist);
            }
            return;
        }

        if (node.right == null) {
            bestChild = node.left;
        } else if (node.left == null) {
            bestChild = node.right;
        } else {
            float pV = getDimension(point, dimension), nV = getDimension(node.obj, dimension);
            if (pV < nV) {
                bestChild = node.left;
            } else {
                bestChild = node.right;
            }
        }

        nearestSearch(point, bestChild, bestNodes, maxNodes);

        if (bestNodes.size() < maxNodes || ownDist < (Float) bestNodes.peek().get("dist")) {
            saveNode(bestNodes, maxNodes, node, ownDist);
        }
        if (bestNodes.size() < maxNodes || Math.abs(linearDist) < (Float) bestNodes.peek().get("dist")) {
            if (bestChild == node.left) {
                otherChild = node.right;
            } else {
                otherChild = node.left;
            }
            if (otherChild != null) {
                nearestSearch(linearPoint, otherChild, bestNodes, maxNodes);
            }
        }

    }

    private void saveNode(BinaryHeap bestNodes, int maxNodes, Node node, float distance) {
        Map<String, Object> tem = new HashMap<String, Object>();
        tem.put("obj", node);
        tem.put("dist", distance);
        bestNodes.push(tem);
        if (bestNodes.size() > maxNodes)
            bestNodes.pop();
    }
    
    public float balanceFactor() {
        return  height(root) / (float) Math.log(count(root)) / (float) Math.log(2);
    }

    private float height(Node node) {
        if (node == null)
            return 0;
        return Math.max(height(node.left), height(node.right)) + 1;
    }
    
    private float count(Node node) {
        if (node == null)
            return 0;
        return count(node.left) + count(node.right) + 1;
    }
}
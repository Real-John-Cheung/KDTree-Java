import kdtree.*;

public class Point{
  public float x, y;
  public Point(float x, float y){
    this.x = x;
    this.y = y;
  }
  public Point(){
    this(0,0);
  }
}

ArrayList<Point> points;
KDTree tree;
void setup(){
  size(800,800);
  tree = new KDTree(new Point[0], null, new String[]{"x", "y"}, this);
  points = new ArrayList<>();
  for(int i = 0; i < 100; i++){
    points.add(new Point(random(width), random(height)));
    tree.insert(points.get(i));
  }
}

void draw(){
  background(255);
  noStroke();
  fill(0);
  Point correct = null;
  float minDist = Float.POSITIVE_INFINITY;
  for(int i = 0; i < points.size(); i++){
    float dist = (points.get(i).x - mouseX) * (points.get(i).x - mouseX) + (points.get(i).y - mouseY) * (points.get(i).y - mouseY);
    if (dist < minDist){
      minDist = dist;
      correct = points.get(i);
    }
    circle(points.get(i).x, points.get(i).y, 10);
  }
  
  fill(0,255,0);
  circle(correct.x, correct.y, 7);
  fill(255,0,0);
  Point tem = new Point(mouseX, mouseY);
  Point nearest = (Point) tree.nearest(tem, 1, -1)[0].get("obj");
  circle(nearest.x, nearest.y, 5);
}

void mouseClicked(){
  if (mouseButton == LEFT){
    Point tem = new Point(mouseX, mouseY);
    Point nearest = (Point) tree.nearest(tem, 1, -1)[0].get("obj");
    tree.remove(nearest);
    points.remove(nearest);
  } else {
    Point tem = new Point(mouseX, mouseY);
    tree.insert(tem);
    points.add(tem);
  }
}

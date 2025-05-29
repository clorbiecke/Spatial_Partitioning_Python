from pygame import draw
import physics_objects, contact
from physics_objects import Bounds, PhysicsObject
from pygame import Vector2
from collections import defaultdict
import math

# QUAD TREE NODE
# Initialized with an anchor position, size, and capacity.
class QuadTreeNode(Bounds):
    def __init__(self, pos:Vector2, size:Vector2, depth:int, max_depth:int, capacity:int=4):
        # print(f"QuadTreeNode... pos: {pos}, size: {size}")
        super().__init__(pos, size)
        # print(f"node depth: {depth}")
        self.capacity = capacity
        self.depth = depth
        self.max_depth = max_depth
        self.divided = False
        self.objects = set()

    # PROPERTIES

    # METHODS
    # Divides node into 4 child nodes and redistributes objects to children.
    def sub_divide(self) -> bool:
        self.divided = True
        # create 4 child nodes
        size = self.size*0.5
        depth = self.depth+1
        self.nw = QuadTreeNode(self.pos, size, depth, self.max_depth, self.capacity)
        self.ne = QuadTreeNode(self.pos+(size.x,0), size, depth, self.max_depth, self.capacity)
        self.sw = QuadTreeNode(self.pos+(0,size.y), size, depth, self.max_depth, self.capacity)
        self.se = QuadTreeNode(self.pos+size, size, depth, self.max_depth, self.capacity)
        # redistribute current objects into child nodes
        for child in self.get_children():
            for obj in self.objects:
                child.insert(obj)
        self.objects = None

    # Insert an object into this node. Returns 'False' if the object is larger than the bounds of this node. 
    # If this node is full, inserts into child node.
    def insert(self, obj:PhysicsObject) -> bool:
        # check if obj is in bounds of node
        if not self.intersects(obj.bounds): return False
        # if not divided, insert into self
        if not self.divided:
            self.objects.add(obj)
            if len(self.objects) > self.capacity and self.depth < self.max_depth: 
                self.sub_divide()
        # if divided, insert into each each child that intersects
        else:
            for child in [self.ne,self.nw,self.se,self.sw]:
                child.insert(obj)

    # Returns a set of all objects in nodes that are within the given boundry range.
    def query(self, bounds:Bounds, found:set = set()) -> set:
        if self.intersects(bounds):
            if not self.divided:
                found.update(self.objects)
            else:
                for child in self.get_children():
                    child.query(bounds, found)
        return found

    # Returns a set of all leaves within the given bounds
    def query_leaves(self, bounds:Bounds, found:set = set()) -> set:
        if self.intersects(bounds):
            if not self.divided:
                found.add(self)
            else:
                for child in self.get_children():
                    child.query_leaves(bounds, found)
        return found
     
    # Returns 'True' if the object is found within this node or its children, else returns 'False'.
    def contains(self, obj) -> bool:
        # check if obj is in self
        if not self.divided:
            return (obj in self.objects)
        # if not, check if obj is in child nodes
        else:
            for child in self.get_children():
                if obj in child: return True
            return False

    # Returns 'True' if this node's boundries intersect the given range.
    def intersects(self, bounds:Bounds):
        # check if boundry intersects range
        return (bounds.pos.x < self.pos.x+self.size.x and bounds.pos.y <self.pos.y+self.size.y 
                and self.pos.x < bounds.pos.x+bounds.size.x and self.pos.y < bounds.pos.y+bounds.size.y)

    # Returns a tuple containing this node's children if this node is a branch. Returns 'False' if this node is not a branch.
    def get_children(self) -> tuple['QuadTreeNode','QuadTreeNode','QuadTreeNode','QuadTreeNode']:
        if self.divided:
            return (self.ne,self.nw,self.se,self.sw)
        else: return False

    # Draw the boundries of this node. If this is a branch node, calls draw() on its children instead.
    def draw(self, surface):
        if not self.divided:
            # print(f"drawing node at depth {self.depth}...")
            draw.rect(surface, (100,100,200), (self.pos, self.size), 1)
        else:
            for child in self.get_children():
                child.draw(surface)


# Will hold the root of the quad tree.
# TODO: have this hold maximum depth and have child nodes reference this one.
class QuadTreeRoot(QuadTreeNode):
    def __init__(self, pos, size, max_depth = 8, node_capacity = 4):
        super().__init__(pos=pos, size=size, depth=0, max_depth=max_depth, capacity=node_capacity)
    



# SPATIAL HASHING TABLE
class SpatialHashTable(Bounds):
    def __init__(self, pos, size, cell_width, cell_height):
        super().__init__(pos, size)
        self.cell_width = cell_width
        self.cell_height = cell_height
        self.grid = defaultdict(set[PhysicsObject])

    def coord_hash(self, x, y):
        return (x // self.cell_width, y // self.cell_height)
    
    def insert(self, obj:PhysicsObject):
        xmin,ymin = self.coord_hash(obj.bounds.pos)
        xmax,ymax = self.coord_hash(obj.bounds.pos+obj.bounds.size)
        for x in range(xmin, xmax+1):
            for y in range(ymin, ymax+1):
                self.grid[(x,y)].add(obj)

    def query(self, bounds:Bounds, found:set[PhysicsObject] = set()):
        xmin,ymin = self.coord_hash(bounds.pos)
        xmax,ymax = self.coord_hash(bounds.pos+bounds.size)
        for x in range(xmin, xmax+1):
            for y in range(ymin, ymax+1):
                found.update(self.grid[(x,y)])
        return found

    def clear(self):
        self.grid.clear()

    def draw(self, surface):
        draw.rect(surface, (200,100,100), (self.pos, self.size), 1)
        x = self.width // self.cell_width
        y = self.height // self.cell_height
        for n in range(len(x+1)):
            for m in range(len(y+1)):
                draw.rect(surface, (200,100,100), (self.pos+(x*self.cell_width, y*self.cell_height), (self.cell_width,self.cell_height)), 1)










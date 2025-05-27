import pygame
from pygame.constants import *
from pygame.math import Vector2

from screeninfo import screeninfo, get_monitors

import colors
from physics_objects import UniformCircle, Wall, PhysicsObject
import contact
from contact import Contact
from spatial_partitioning import QuadTreeRoot as QuadTree

import math, random, time



# INITIALIZE PYGAME AND OPEN WINDOW
pygame.init()
pygame.font.init()

monitor = get_monitors()[0]
win_w, win_h = (0.5*monitor.width, 0.9*monitor.height)
window = pygame.display.set_mode([win_w, win_h])

# SET UP TIMING
FPS = 60
dt = 1/FPS
clock = pygame.time.Clock()

# COLORS
BG_COLOR = colors.darken(colors.midnight, 0.5)


# QUAD TREE
tree_pos =Vector2(10,10)
tree_size = Vector2(win_w-20, win_h-20)
tree_depth = 6
tree_cap = 4
tree = QuadTree(tree_pos, tree_size)

# OBJECTS
objects:set[PhysicsObject] = set()
wall_pts = [tree_pos, tree_pos+(tree_size.x,0), tree_pos+tree_size, tree_pos+(0,tree_size.y)]
for i in range(len(wall_pts)):
    objects.add(Wall(wall_pts[i-1], wall_pts[i], colors.red, 0, f"BoundryWall{i}", 0))

# Collisions
restitution = 1
friction = 0.5

## Spawning
num_to_spawn = 500
spawn_count = 0
def spawn_circle():
    global spawn_count
    spawn_count += 1
    pos = (random.uniform(tree.pos.x, tree.pos.x+tree.width), random.uniform(tree.pos.y, tree.pos.y+tree.height))
    vel = (random.uniform(-100,100), random.uniform(-100,100))
    objects.add(UniformCircle(density=0.001, radius = 10, pos=pos, vel=vel))



# SETUP SIMULATION
running = True
start = 0
end = 0
def perf_time() -> float:
    return end-start

# SIMULATION LOOP
while running:
    print("START OF LOOP")
    if spawn_count < num_to_spawn:
        for i in range(10):
            spawn_circle()
    print(f"spawn count: {spawn_count}")

    # EVENTS
    while event := pygame.event.poll():
        if event.type == QUIT or (event.type == KEYDOWN and event.key == K_ESCAPE):
            running = False
    
    # PHYSICS UPDATE
    start = time.perf_counter()
    for obj in objects: 
        obj.clear_force()
    
    for obj in objects:
        obj.update(dt)
    end = time.perf_counter()
    print(f"updating objects... {perf_time()}")

    # COLLISION UPDATE
    start = time.perf_counter()
    tree = QuadTree(tree_pos,tree_size, tree_depth, tree_cap)
    for obj in objects:
        tree.insert(obj)
    end = time.perf_counter()
    print(f"creating quad tree... {perf_time()}")

    start = time.perf_counter()
    for obj in objects:
        others = tree.query(obj.bounds)
        others.discard(obj)
        for other in others:
            c:Contact = contact.generate(obj, other, restitution=restitution, friction=friction)
    end = time.perf_counter()
    print(f"resolving collisions... {perf_time()}")

    # DRAW ALL
    start = time.perf_counter()
    window.fill(BG_COLOR)
    for obj in objects:
        obj.draw(window)
    tree.draw(window)
    end = time.perf_counter()
    print(f"drawing objects... {perf_time()}")

    # UPDATE DISPLAY
    pygame.display.update()
    clock.tick(FPS)
    print("")

    if spawn_count >= num_to_spawn: running = False

# END OF SIMULATION LOOP
pygame.quit()
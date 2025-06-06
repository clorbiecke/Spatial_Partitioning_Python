import pygame
from pygame.constants import *
from pygame.math import Vector2

from screeninfo import screeninfo, get_monitors

import colors
from physics_objects import UniformCircle, Wall, PhysicsObject
import contact
from contact import Contact
from spatial_partitioning import *

import math, random, time, itertools




# INITIALIZE PYGAME AND OPEN WINDOW
pygame.init()
pygame.font.init()

monitor = get_monitors()[0]
win_w, win_h = (0.9*monitor.width, 0.9*monitor.height)
window = pygame.display.set_mode([win_w, win_h])

# SET UP TIMING
FPS = 60
dt = 1/FPS
clock = pygame.time.Clock()

# COLORS
BG_COLOR = colors.darken(colors.midnight, 0.5)


# QUAD TREE
tree_pos = Vector2(10,10)
tree_size = Vector2(win_w-20, win_h-20)
tree_depth = 8
tree_cap = 6
tree = QuadTreeRoot(tree_pos, tree_size)

# SPATIAL HASH TABLE
grid_pos = Vector2(10,10)
grid_size = Vector2(win_w-20, win_h-20)
cell_width = grid_size.x / 10
cell_height = grid_size.y / 10
# hash_grid = SpatialHashTable(grid_pos, grid_size, cell_width, cell_height)

# OBJECTS
objects:set[PhysicsObject] = set()
# objects:list[PhysicsObject] = []

wall_pts = [tree_pos, tree_pos+(tree_size.x,0), tree_pos+tree_size, tree_pos+(0,tree_size.y)]
for i in range(len(wall_pts)):
    objects.add(Wall(wall_pts[i-1], wall_pts[i], colors.red, 0, f"BoundryWall{i}", 0))
    # objects.append(Wall(wall_pts[i-1], wall_pts[i], colors.red, 0, f"BoundryWall{i}", 0))


# Collisions
restitution = 1
friction = 0.5

## Spawning
num_to_spawn = 10000
spawn_count = 0
def spawn_circle():
    global spawn_count
    spawn_count += 1
    pos = (random.uniform(tree.pos.x, tree.pos.x+tree.width), random.uniform(tree.pos.y, tree.pos.y+tree.height))
    vel = (random.uniform(-100,100), random.uniform(-100,100))
    objects.add(UniformCircle(density=0.001, radius = 0.5, pos=pos, vel=vel))
    # objects.append(UniformCircle(density=0.001, radius = random.randrange(1,20), pos=pos, vel=vel))



# SETUP SIMULATION
running = True
start = 0
end = 0


def setup_sim(num_iterations:int = 1):
    global running, start, end, sim_num, sim_iter, run_times
    running = True
    start = 0
    end = 0
    sim_num = 0
    sim_iter = num_iterations
    run_times = []

    reset_sim()

def reset_sim():
    # set/reset objects and spawning
    global spawn_count, num_to_spawn
    objects.clear()
    spawn_count = 0

def next_sim() -> bool:
    global sim_end, sim_start, sim_num
    sim_end = time.perf_counter()
    sim_num += 1
    end_sim()
    reset_sim()
    sim_start = time.perf_counter()
    return sim_num < sim_iter

def start_sim():
    global sim_start
    # setup_sim()
    sim_start = time.perf_counter()

def end_sim():
    run_times.append(sim_end-sim_start)
    # print(f"* Sim iteration {sim_num}\{sim_iter} ended. Time to complete: {sim_end-sim_start}\n")

def get_avg_time():
    avg = 0
    for time in run_times:
        avg += time
    return (avg/len(run_times))


def perf_time() -> float:
    return end-start

# UPDATE SPATIAL PARTITION
def update_tree():
    global tree
    tree = QuadTreeRoot(tree_pos,tree_size, tree_depth, tree_cap)
    # print(f"updating... num leaves: {len(tree.query_leaves(tree))}, should be 1")

    for obj in objects:
        tree.insert(obj)
    
def update_hash_grid():
    global hash_grid
    hash_grid.clear()
    for obj in objects:
        hash_grid.insert(obj)

# PROCESS COLLISIONS
# quad tree
def process_all_pairs():
    global start, end
    num_contacts = 0
    start = time.perf_counter()
    for a,b in itertools.combinations(objects, 2):
        c:Contact = contact.generate(a, b)
        if c.bool:
            num_contacts += 1
            c.resolve(restitution=restitution, friction=friction)
    end = time.perf_counter()
    print(f"process_all_pairs() | num_contacts: {num_contacts}, time: {perf_time()}")

def process_quadtree():
    global start, end
    num_contacts = 0
    start = time.perf_counter()
    for obj in objects:
        others = tree.query(obj.bounds)
        others.discard(obj)
        for other in others:
            c:Contact = contact.generate(obj, other)
            if c.bool:
                num_contacts += 1
                c.resolve(restitution=restitution, friction=friction)
    end = time.perf_counter()
    print(f"process_quadtree() | num_contacts: {num_contacts}, time: {perf_time()}")

def process_quadtree_no_duplicates():
    global start, end
    num_contacts = 0
    start = time.perf_counter()
    checked = set()
    for obj in objects:
        checked.add(obj)
        others = tree.query(obj.bounds)
        for other in others:
            if other in checked: continue
            c:Contact = contact.generate(obj, other)
            if c.bool:
                num_contacts += 1
                c.resolve(restitution=restitution, friction=friction)
    end = time.perf_counter()
    print(f"process_quadtree_no_duplicates() | num_contacts: {num_contacts}, time: {perf_time()}")

def process_quadtree_leaves():
    # get all leaves in quad tree and check collisions in each
    global start, end
    num_contacts = 0
    num_checks = 0
    start = time.perf_counter()
    leaves:set[QuadTreeNode] = tree.query_leaves(tree)
    # print(f"leaves len: {len(leaves)}")
    for leaf in leaves:
        for a,b in itertools.combinations(leaf.objects, 2):
            num_checks += 1
            c:Contact = contact.generate(a,b)
            if c.bool:
                num_contacts += 1
                c.resolve(restitution=restitution, friction=friction)
    end = time.perf_counter()
    print(f"process_quadtree_leaves() | num_contacts: {num_contacts}, time: {perf_time()}")#, num_checks: {num_checks}"

# hash table
# def process_hash_table():
#     global start, end
#     num_contacts = 0
#     start = time.perf_counter()
#     for obj in objects:
#         found = hash_grid.query(obj.bounds)
#         found.discard(obj)
#         for other in found:
#             c:Contact = contact.generate(obj, other)
#             if c.bool:
#                 num_contacts += 1
#                 c.resolve(restitution=restitution, friction=friction)
#     end = time.perf_counter()
#     print(f"process_hash_table() | num_contacts: {num_contacts}, time: {perf_time()}")   

# def process_hash_table_cells():
#     global start, end
#     num_contacts = 0
#     start = time.perf_counter()
#     for cell in hash_grid:
#         for a,b in itertools.combinations(cell, 2):
#             c:Contact = contact.generate(a,b)
#             if c.bool:
#                 num_contacts += 1
#                 c.resolve(restitution=restitution, friction=friction)
#     end = time.perf_counter()
#     print(f"process_hash_table_cells() | num_contacts: {num_contacts}, time: {perf_time()}")   


# DEBUGGING
setup_sim(1)
processing_times = []
# END DEBUGGING

# SIMULATION LOOP
start_sim()
while running:
    print("START OF LOOP")
    if spawn_count < num_to_spawn:
        for i in range(num_to_spawn // 100):
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
    ## update quadtree
    start = time.perf_counter()
    update_tree()
    end = time.perf_counter()
    print(f"creating quad tree... {perf_time()}")

    ## update hash table
    # start = time.perf_counter()
    # update_hash_grid()
    # end = time.perf_counter()
    # print(f"creating hash table... {perf_time()}")

    ## process collisions
    # start = time.perf_counter()
    # process_all_pairs()
    # process_quadtree()
    # process_quadtree_no_duplicates()
    process_quadtree_leaves()
    processing_times.append(perf_time())
    # process_hash_table()
    # process_hash_table_cells()
    # end = time.perf_counter()
    # print(f"resolving collisions... {perf_time()}")

    # DRAW ALL
    start = time.perf_counter()
    window.fill(BG_COLOR)
    for obj in objects:
        obj.draw(window)
    tree.draw(window)
    # hash_grid.draw(window)
    end = time.perf_counter()
    print(f"drawing objects... {perf_time()}")

    # UPDATE DISPLAY
    pygame.display.update()
    clock.tick(FPS)
    # print("")

    if spawn_count >= num_to_spawn: 
        running = next_sim()

# END OF SIMULATION LOOP
print(f"\nSimulations complete. Average iteration time: {get_avg_time()}")
avg_proc_time = 0
for t in processing_times:
    avg_proc_time += t
avg_proc_time /= len(processing_times)
print(f"Average processing time (capacity: {tree_cap}): {avg_proc_time}\n")
pygame.quit()
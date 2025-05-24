import pygame
from pygame.constants import *
from pygame.math import Vector2, Vector3
from physics_objects import Ball
from screeninfo import get_monitors
from forces import Gravity, AirDrag
from hud_components import Bar
import math
import random
#from physics_objects import Circle  # copy in physics_objects.py from your previous project
from forces import *

# INITIALIZE PYGAME AND OPEN WINDOW
pygame.init()
monitor = get_monitors()[0]
screen_width, screen_height = (monitor.width, monitor.height)
window = pygame.display.set_mode([0.8*screen_width, 0.8*screen_height])
print(f"Window size: {window.get_size()}")

# FONTS
font = pygame.font.SysFont('impact', 20)
wind_font = pygame.font.SysFont('impact', int(min(window.get_size())/12))

# SCALE
PPM = screen_height/4
MPP = 1.0/PPM
print(f"PPM: {PPM}")

# SETUP TIMING
fps = 60*3
dt = 1/(fps*1)
clock = pygame.time.Clock()

# COLORS
BLACK = (0,0,0,128)
SKY_BLUE = (100,100,255,255)
WHITE = (255,255,255,128)



# SETUP OBJECTS
objects = [] 
## Circles
num_balls = 15
ball_radius = 0.035 * PPM
ball_spacing = 0.105 * PPM
ball_mass = .010 #kilograms
# Ball.DRAG_C *= MPP**2
pairs = []

BALL_COLS = [(255,0,0), (0,255,0), (255,180,0)]
cols_it = itertools.cycle(BALL_COLS)
# step_col = int((255 - 55)/(num_balls/len(BALL_COLS)))
num_steps = int(num_balls/len(BALL_COLS))+1

for i in range(num_balls):
    ncol = next(cols_it)
    k = (int(i/len(BALL_COLS)))/num_steps
    col = tuple(c - int(k * (3*c/5)) for c in ncol)
    ball = Ball(mass=ball_mass, pos=(window.get_width()/2, ball_spacing * (i + 1)), radius=ball_radius, color=col, static=(i == num_balls-1), name=f"Ball {i}")  # change mass and radius
    objects.append(ball)
    if i > 0: pairs.append((objects[len(objects)-2], objects[len(objects)-1]))

# SETUP FORCES
## gravity
gravity = Gravity(objects=objects, gravity=(0,9.8 * PPM))
## spring
stiffness = 10000*MPP
damping = (2*math.sqrt(stiffness*ball_mass))*0.5
spring = SpringForce(length=ball_spacing, stiffness=stiffness, damping=damping, pairs=pairs)
## wind and drag
drag = AirDrag(air_density=1.225*(MPP**3), wind_vel=(0,0), objects=objects)
## repulsion when balls touch (makes sense, same)
repulsion = SpringRepulsion(10, objects=objects)

# particles
particles = []
max_particles = 10000
ppf = 100
offset = 20
min_size, max_size = (4, 12)
lifespan = 1000
delta = 0.01
par_color = (128,190,128,128)
par_color_to = (-128,-190,-128,-128)
def par_vel():
    vel = Vector2(0)  
    vel.x = drag.wind_vel.x
    vel.y = random.randrange(-int(abs(drag.wind_vel.x/2)), int(abs(drag.wind_vel.x))) if abs(drag.wind_vel.x) >= 1 else 0
    return vel
def par_delta(particle):
    global delta, lifespan
    if particle.size < min_size:
        particle.curr_delta = delta
        return particle.curr_delta
    elif particle.size > max_size: 
        particle.curr_delta = -delta
        return particle.curr_delta
    else: return particle.curr_delta
def spawn_particles():
    global particles, ppf, offset, drag, window, dt
    if len(Particle.particles) > max_particles: return
    for i in range(ppf):
        pos = Vector2(0)
        pos.x = random.randrange(-offset,int(window.get_width()+offset))
        pos.y = random.randrange(-offset,int(window.get_height()+offset))
        vel = Vector2(0)  
        vel.x = drag.wind_vel.x
        vel.y = random.randrange(-int(drag.wind_vel.magnitude()/2), int(drag.wind_vel.magnitude())) if abs(drag.wind_vel.x) >= 1 else 0
        # p = Particle(random.randrange(80,120), par_color, 3000, 0, par_color_to, pos=pos, vel=lambda: par_vel(), mass=0.001, delta_size=lambda: par_delta(p))
        p = Particle(random.randrange(min_size,max_size), par_color, lifespan, 0, par_color_to, pos=pos, vel=lambda: par_vel(), mass=0.001)
        # p.curr_delta = delta
        # p.delta_size=lambda: par_delta(p)
def draw_particles(window):
    Particle.draw_all(window)
def draw_particle_txt(window):
    if Particle.MAX_PARTICLES != 0:
        txt = font.render("Press '0' to stop particle spawning.", True, WHITE)
        window.blit(txt, (window.get_width()-txt.get_width()-20, 10))

# HOLDING BALLS
held = None
def ball_at_mouse(pos):
    for obj in objects:
        if (obj.pos - pos).magnitude() <= ball_radius:
            return obj
    return None
def delete_held():
    global objects, pairs, held
    objects.remove(held)
    n = len(pairs)
    for i in range(n):
        if pairs[n-(i+1)].__contains__(held):
            pairs.pop(n-(i+1))

    held = None

# DEBUGGING
show_debug = False
debug_color = (255,255,255)

def toggle_debug():
    show_debug = not show_debug
def init_debug():
    debug = []
    debug.append(f"Held object: {held.name if held is not None else "None"}")
    debug.append(f"Wind Velx: {drag.wind_vel.x}")
    debug.append(f"curr/max: {wind_bar.curr_value}/{wind_bar.max_value}")
    debug.append(f"fill pos/size: {wind_bar.fill_rect}")
    debug.append(f"fullness: {wind_bar.fullness}")
    return debug
def draw_debug(window):
    debug = init_debug()
    for i in range(len(debug)):
        txt = font.render(debug[i], True, debug_color)
        window.blit(txt, (0, i * txt.get_height()))
# HUD
bar_width, bar_height = (window.get_width()/3, window.get_height()/10)
wind_bar = Bar((window.get_width()/2 - bar_width/2, 10), (bar_width, bar_height), AirDrag.MAX_WIND_SPD, drag.wind_vel.x, 'center', 'horizontal', True, BLACK, 0, WHITE)
wind_bar.border_width=0
wind_txt = ""
def draw_hud(window):
    wind_bar.set_curr_value(drag.wind_vel.x)
    wind_txt = f"{round(drag.wind_vel.x/PPM, 1)}m/s"
    wind_render = wind_font.render(wind_txt, True, WHITE)
    wind_render_pos = Vector2(wind_bar.pos.x-(5+wind_render.get_width()), wind_bar.pos.y + (wind_bar.size.y-wind_render.get_height())/2)
    wind_bar.draw(window)
    window.blit(wind_render, wind_render_pos)
    draw_particle_txt(window)
    if show_debug:
        draw_debug(window)

# STATES

## pause state
def toggle_pause_state():
    global state
    if state == "running":
        state = "paused"
        for obj in objects:
            obj.vel = Vector2(0)
    else: state = "running"
## add mode
in_add_mode = False
num_added = 0
prev_ball = None
def add_mode_off():
    global num_added, prev_ball, in_add_mode
    num_added = 0
    prev_ball = None
    in_add_mode = False
def add_mode_on():
    global in_add_mode
    print("Add mode on!")
    in_add_mode = True
def add_ball(pos):
    global objects, num_added, prev_ball, pairs
    col=next(cols_it)
    static=(prev_ball is None)
    name=f"New Ball{num_added}"
    print(f"prev_ball: {prev_ball.name if prev_ball is not None else "None"}")
    ball = Ball(mass=ball_mass, pos=pos, radius=ball_radius, color=col, static=static, name=name)
    objects.append(ball)
    num_added += 1
    if prev_ball is not None:
        pairs.append((prev_ball, ball))
    prev_ball = ball
def link_to_prev(ball):
    global pairs, prev_ball
    if prev_ball is None: return
    for p in pairs:
        if p.__contains__(ball) and p.__contains__(prev_ball):
            print(f"A bond already exists between these balls. (lol)")
            return
    pairs.append((prev_ball, ball))
def del_link_to_prev(ball):
    global prev_ball, pairs
    n = len(pairs)
    for i in range(n):
        if pairs[i].__contains__(ball) and pairs[i].__contains__(prev_ball):
            pairs.pop(i)
            return
def draw_selected(window):
    if prev_ball is not None:
        pygame.draw.circle(window, prev_ball.color, prev_ball.pos, prev_ball.radius + 3, width=1)
        pygame.draw.circle(window, WHITE, prev_ball.pos, prev_ball.radius + 2, width=1)
# game loop
state = "running"
clock.tick()
while state != "quit":
    # DISPLAY
    pygame.display.update()
    # TIMING
    clock.tick(fps)
    # BACKGROUND GRAPHICS
    window.fill(SKY_BLUE)
    spawn_particles()
    # EVENTS
    while event := pygame.event.poll():
        if event.type == QUIT or event.type == KEYDOWN and event.key == K_ESCAPE:
            state = "quit"
        if event.type == KEYDOWN:
            if event.key == K_SPACE and held is not None:
                held.toggle_static()
            elif event.key == K_p:
                toggle_pause_state()
            elif event.key == K_F4:
                toggle_debug()
            elif event.key == K_LSHIFT:
                add_mode_on()
            elif event.key == K_0:
                Particle.MAX_PARTICLES = 0
                max_particles = 0
                ppf = 0
        if event.type == KEYUP:
            if event.key == K_LSHIFT:
                add_mode_off()
        if in_add_mode:
            if event.type == MOUSEBUTTONDOWN and event.button == 1:
                # if no balls have been added yet
                clicked = ball_at_mouse(event.pos)
                if prev_ball is None:
                    prev_ball = clicked
                    if prev_ball is None:
                        add_ball(event.pos)
                # if a ball has already been added
                elif clicked is not prev_ball:
                    if clicked is None:
                        add_ball(event.pos)
                    else:
                        link_to_prev(clicked)
                        prev_ball = clicked
            elif event.type == MOUSEBUTTONDOWN and event.button == 3:
                clicked = ball_at_mouse(event.pos)
                if clicked is not None:
                    if prev_ball is not None:
                        print("Calling del_link_to_prev(clicked)...")
                        del_link_to_prev(clicked)
    
    if state == "running":
        # PHYSICS
        ## clear all forces from each object
        for o in objects:
            o.clear_force()
        for p in particles:
            p.clear_force()

        ## apply all forces
        gravity.apply()
        spring.apply()
        drag.apply()
        repulsion.apply()


        ## Change the following to act on all objects using a for loop
        ## update all objects
        for o in objects:
            o.update(dt)
        Particle.update_all(dt)    
    # STATE CHECKS
    ## Mouse state check for grabbing objects
    mouse = pygame.mouse.get_pressed()
    key = pygame.key.get_pressed()
    if mouse[0]: 
        mpos = Vector2(pygame.mouse.get_pos())
        if held is None:
            held = ball_at_mouse(mpos)
        elif held is not None:
            # There's no point in keeping track of the prev mouse pos, because it will be somewhere within radius, so just use ball pos
            # The ball being at a constant offset is gross
            if state == "running":
                held.vel =  (mpos - held.pos)/dt
            elif state == "paused":
                held.pos = mpos
        if key[K_LCTRL] and held is not None:
            delete_held()
    elif held is not None:
        if held.static:
            held.vel = Vector2(0)
        held = None
    ### Add Mode: 

    ## Key state check for changing wind velocity
    if key[K_LEFT]:
        drag.dec_wind_velx()
        for p in particles:
            p.vel.x = drag.wind_vel.x
    if key[K_RIGHT]:
        drag.inc_wind_velx()
        for p in particles:
            p.vel.x = drag.wind_vel.x

    # GRAPHICS
    ## draw all objects
    for o in objects:
        o.draw(window)
    spring.draw(window)
    draw_particles(window)
    if in_add_mode: draw_selected(window)
    draw_hud(window)


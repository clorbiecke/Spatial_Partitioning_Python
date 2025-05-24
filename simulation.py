import pygame
from pygame.constants import *
from pygame.math import Vector2

from screeninfo import screeninfo, get_monitors

import colors

import math, random



# INITIALIZE PYGAME AND OPEN WINDOW
pygame.init()
pygame.font.init()

monitor = get_monitors()[0]
win_w, win_h = (0.9*monitor.width, 0.9*monitor.height)
window = pygame.display.set_mode([win_w, win_h])

# SET UP TIMING
FPS = 60
clock = pygame.time.Clock()

# COLORS
BG_COLOR = colors.darken(colors.midnight, 0.5)










# SETUP SIMULATION
running = True

# SIMULATION LOOP
while running:
    # EVENTS
    while event := pygame.event.poll():
        if event.type == QUIT or (event.type == KEYDOWN and event.key == K_ESCAPE):
            running = False

    # DRAW ALL
    window.fill(BG_COLOR)

    # UPDATE DISPLAY
    pygame.display.update()
    clock.tick(FPS)

# END OF SIMULATION LOOP
pygame.quit()
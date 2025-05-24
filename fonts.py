import pygame
from pygame.locals import *
import pygame.font as font

pygame.font.init()

txt_20 = font.SysFont('roboto', 20)
txt_30 = font.SysFont('roboto', 30)
txt_40 = font.SysFont('roboto', 40)
txt_50 = font.SysFont('roboto', 50)
txt_60 = font.SysFont('roboto', 60)
def txt_size(size:int, font_name:str='roboto'):
    return font.SysFont(font_name, int(size))

head_60 = font.SysFont('impact', 60, bold=False)
head_80 = font.SysFont('impact', 80, bold=False)
head_100 = font.SysFont('impact', 100, bold=False)
head_120 = font.SysFont('impact', 120, bold=False)
def head_size(size:int, font_name:str='impact', bold=False):
    return font.SysFont(font_name, int(size), bold)
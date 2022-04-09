#seed_gen.py
from random import seed
from random import randint

upp = 2**31-1
until = int(1e4)

seed(19890421)
for i in range(until):
  print(randint(0,upp))

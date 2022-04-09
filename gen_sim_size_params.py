# gen_sim_size_params.py

resp = 40
non_resp = 20

for i in range(1,non_resp+1):
  for j in range(1,resp+1):
    if i+j >= 8:
      print(i,j)

setwd('G:/博士/泛基因组组/PGG-project/analysis/SV stats/SV_ld/02_LD_common_rare/')

read.table("all.LD_all_common")->all_common;
read.table("all.LD_all_rare")->all_rare;
p1=ggplot() +
  geom_line(data = all_common, aes(x = all_common[,1]/1000, y = all_common[,2]), color = "black", size = 1) +
  geom_line(data = all_rare, aes(x =  all_rare[,1]/1000, y = all_rare[,2]), color = "black", size = 1, linetype = "dashed")+labs(
    x = "Distance (Kb)",
    y = expression("mean LD (r"^2*")"),
    title='ALL SVs'
  )+theme_classic()+scale_x_continuous(limits = c(0, 100))

read.table("del.LD_del_common")->del_c;
read.table("del.LD_del_rare")->del_r;
p2=ggplot() +
  geom_line(data = del_c, aes(x = del_c[,1]/1000, y = del_c[,2]), color = "black", size = 1) +
  geom_line(data = del_r, aes(x =  del_r[,1]/1000, y = del_r[,2]), color = "black", size = 1, linetype = "dashed")+labs(
    x = "Distance (Kb)",
    y = expression("mean LD (r"^2*")"),
    title='DEL'
  )+theme_classic()+scale_x_continuous(limits = c(0, 100))

read.table("ins.LD_INS_common")->ins_l;
read.table("ins.LD_INS_rare")->ins_r;
p3=ggplot() +
  geom_line(data = ins_l, aes(x = ins_l[,1]/1000, y = ins_l[,2]), color = "black", size = 1) +
  geom_line(data = ins_r, aes(x =  ins_r[,1]/1000, y = ins_r[,2]), color = "black", size = 1, linetype = "dashed")+labs(
    x = "Distance (Kb)",
    y = expression("mean LD (r"^2*")"),
    title='INS'
  )+theme_classic()+scale_x_continuous(limits = c(0, 100))


read.table("dup.LD_dup_common")->dup_c;
read.table("dup.LD_dup_rare")->dup_r;
p4=ggplot() +
  geom_line(data = dup_c, aes(x = dup_c[,1]/1000, y = dup_c[,2]), color = "black", size = 1) +
  geom_line(data = dup_r, aes(x =  dup_r[,1]/1000, y = dup_r[,2]), color = "black", size = 1, linetype = "dashed")+labs(
    x = "Distance (Kb)",
    y = expression("mean LD (r"^2*")"),
    title='DUP'
  )+theme_classic()+scale_x_continuous(limits = c(0, 100))

read.table("inv.LD_INV_common")->inv_c;
read.table("inv.LD_INV_rare")->inv_r;
p5=ggplot() +
  geom_line(data = inv_c, aes(x = inv_c[,1]/1000, y = inv_c[,2]), color = "black", size = 1) +
  geom_line(data = inv_r, aes(x =  inv_r[,1]/1000, y = inv_r[,2]), color = "black", size = 1, linetype = "dashed")+labs(
    x = "Distance (Kb)",
    y = expression("mean LD (r"^2*")"),
    title='INV'
  )+theme_classic()+scale_x_continuous(limits = c(0, 100))+scale_y_continuous(limits = c(0, 1))

read.table("bnd.LD_bnd_common")->bnd_c;
read.table("bnd.LD_bnd_rare")->bnd_r;
p6=ggplot() +
  geom_line(data = bnd_c, aes(x = bnd_c[,1]/1000, y = bnd_c[,2]), color = "black", size = 1) +
  geom_line(data = bnd_r, aes(x =  bnd_r[,1]/1000, y = bnd_r[,2]), color = "black", size = 1, linetype = "dashed")+labs(
    x = "Distance (Kb)",
    y = expression("mean LD (r"^2*")"),
    title='BND'
  )+theme_classic()+scale_x_continuous(limits = c(0, 100))+scale_y_continuous(limits = c(0, 1))

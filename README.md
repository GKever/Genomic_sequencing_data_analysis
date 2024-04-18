这里主要记录一些我用于基因组测序数据分析的一些笔记和代码，有空的话我会不断地进行一些更新。

一、搭建分析平台
![20240325-1](https://github.com/GKever/Genomic_sequencing_data_analysis/assets/111635048/6bbf124c-aca9-4e18-b781-235fa99e85d9)

测序分析过程基本都时基于Linux平台，可以在服务器上部署也可以在windows/Mac电脑上部署，强烈建议在Windows中使用WSL功能，安装Linux子系统，本文后续都是基于WSL平台（ubuntu 18.04)进行。具体方法为：
1. windows启用wsl功能
   win11中在设置中，系统->可选功能->更多windows功能->勾选“适用于linux的windows子系统”,“虚拟机平台”，“Windows虚拟机监控程序平台”，确定后重启电脑。
2. 安装Linux分发版本
   在Windows商店中搜索Linux，安装自己喜欢的Linux版本，我个人常用Ubuntu系统。
3. Linux系统中搭建分析流程
   启动Linux系统，安装conda软件（https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html），
mkdir -p ~/miniconda3
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh
bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
rm -rf ~/miniconda3/miniconda.sh

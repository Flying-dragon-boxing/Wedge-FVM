**English Version**: 

This is a finite volume method solver for supersonic flow over a wedge based on structured mesh. Van Leer, Roe and KFVS scheme have been used for calculating the fluxes. 

# 计算流体力学 期末大作业 程序文档

## 1. 项目依赖

程序的网格生成和后处理部分依赖于Python3的numpy和matplotlib库, 请确保安装了Python3和这两个库. 我们的程序在Python 3.10.14下测试通过.

程序的计算部分由C++语言写成, 且依赖于xmake, 请确保运行环境中有可用的C++编译器和xmake. 我们的程序在Windows 11 (Visual Studio 2022)和Archlinux (g++)下测试通过.

## 2. 项目结构

C++部分: 

- `main.cpp`: 主程序入口
- `Mesh.h`: 主要的头文件
- `Mesh.cpp`: 处理网格的源文件
- `Flux.cpp`: 计算通量的源文件

Python部分:

- `MeshGen.py`: 生成网格的Python脚本, 在帮助文档给出的网格划分程序的基础上修改而来
- `Plot.py`: 后处理的Python脚本

## 3. 使用方法

### 3.1 网格生成

在`MeshGen.py`中修改网格生成的参数, 运行

```bash
python MeshGen.py
```

即可生成网格文件. 生成的网格文件在`Mesh`文件夹下.

### 3.2 编译

```bash
xmake
```

当编译器路径在环境变量中时, 这应当可以直接解决依赖问题并编译出可执行文件.

### 3.3 计算

```bash
xmake run HW6 <迭代次数> <通量处理方法>
```

其中`<迭代次数>`是一个整数, `<通量处理方法>`是一个字符串, 表示通量的处理方法, 可以是`roe`(Roe格式), `vl`(Van Leer格式)或`kfvs`(KFVS格式).

计算结果会以Numpy的二进制格式保存在项目根目录中. 

### 3.4 后处理

```bash
python Plot.py
```

绘制出计算结果的图像, 并保存在项目根目录中. 此外, 后处理脚本还会输出计算结果的压强比和马赫数.

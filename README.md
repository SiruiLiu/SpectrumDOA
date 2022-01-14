# 算法说明

## 1.概述

    该算法程序主要包含MUSIC核心测向算法，通道校正算法及配置文件读写程序。

## 2.文件内容

    压缩包中主要包含Config文件夹，src文件夹，SimData及CMakeLists.txt.

| 文件或文件夹         | 说明      |
|:--------------:|:-------:|
| Config         | 包含配置文件  |
| src            | 头文件及源文件 |
| SimData        | 仿真调试用数据 |
| CMakeLists.txt | CMake文件 |

在src文件夹中包含另外三个文件夹及mian.cpp文件，说明如下

| 文件夹              | 说明        |
|:----------------:|:---------:|
| MUSIC_Algorithm  | MUSIC测向算法 |
| PhaseCalibration | 通道校正算法    |
| Configure        | 配置文件读写算法  |
| main.cpp         | 仿真测试用主程序  |

## 3.配置文件说明

    配置文件定义了该部分算法模块所需的参数及通道校正参数，具体说明如下

##### 1.AntennaNum

使用天线个数，亦即通道个数，最大不超过8；

##### 2.ArrayRadius

天线阵半径。因天线阵为圆阵，故通过给定的半径和最大通道个数可以确定各天线位置；

##### 3.DataLength

数据长度，用于在初始化时确定用于保存数据的内存空间的大小，且在信源个数估计中及部分其他计算中也需要该参数进行计算；

##### 4.AzimuthStep

方位角搜索精度，单位为度（°）。

##### 5.PitchStep

俯仰角搜索精度，单位为度（°）。

两个方位上搜索精度越高（值越小），运算时间越长。

##### 6.START_FREQ

扫频起始频率

##### 7.STOP_FREQ

扫频截止频率

##### 8.FREQSTEP

扫频步进，以上三个参数定义了频表，即进行了通道校正的频率

##### 9.校正表

后续键名为频率，键值为校正参数。

## 4.使用说明

Configure类中主要为算法相关参数，对应的配置文件读写也在算法程序中进行，因此只需要将头文件包含到主程序中即可；

主要使用类为PhaseCalibration和MUSIC两个类，说明如下：

#### 1.PhaseCalibration

```cpp
/*
 *声明类，参数DataAddr为指向std::vector<complex<float>>类型的指针
 */
PhaseCalibration Calibrator(DataAddr);
/*
 *根据对应频点的采集数据进行校正，Fc为对应频点，该函数会将对应频点及
 *计算所得校正数据写入Config.ini文件中
 */
Calibrator.Calibration(Fc);
```

#### 2.MUSIC

```cpp
/*
 *声明类，参数DataAddr为指向std::vector<complex<float>>类型的指针，
 *参数Fc为AD采样的中心频率
 */
MUSIC Music(DataAddr, Fc);
/*
 *运行算法,Azimuth为计算所得俯仰角序列，Pitch为计算所得方位角序列，两
 *序列长度相同，序列长度即为单频信号个数。
 */
Music.run(Azimuth, Pitch);
```

具体使用可以还可以参考main.cpp程序，程序从SimData文件夹中读取8通道数据，并存入预先声明的内存空间中，并进行测向及及校正等操作。

**<mark>要运行main.cpp，Config.ini配置文件默认配置如下</mark>**

```ini
[AntennaArray]
AntennaNum=8
ArrayRadius=100
[DATA]
DataLength=4096
[SEARCH_SETTING]
AzimuthStep=0.5
PitchStep=0.5
[CALI_TABLE]
START_FREQ=1000
STOP_FREQ=3000
FREQSTEP=100
1200=2.097913,0.343672,2.724445,-0.210563,-2.308476,-0.554236,-2.935008
1300=2.097913,0.343672,2.724445,-0.210563,-2.308476,-0.554236,-2.935008
```

<mark>另，最终数据类型需要根据前端数据类型决定，需要进行优化。</mark>



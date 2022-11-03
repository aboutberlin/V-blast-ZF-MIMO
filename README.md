# V-blast ZF & MIMO
Implemented QR，SORTED QR，MMSE，MMSE QR，V-BLAST，ML algorithm

Implemented a matrix class, through operator overloading for input and output of matrices, addition, subtraction, multiplication and division,slicing, QR decomposition, improved Gram-Schmidt orthogonality and matrix inversion, etc.

This cpp is used to implement paper 

1：Wübben, Dirk & Rinas, Jürgen & Boehnke, R. & Kuehn, Volker & Kammeyer, K.. (2002). Efficient Algorithm for Detecting Layered Space-Time Codes. http://www.ant.uni-bremen.de/pub/papers/wrbkk02_paper.pdf. 

2：D. Wubben, R. Bohnke, V. Kuhn and K. . -D. Kammeyer, "MMSE extension of V-BLAST based on sorted QR decomposition," 2003 IEEE 58th Vehicular Technology Conference. VTC 2003-Fall (IEEE Cat. No.03CH37484), 2003, pp. 508-512 Vol.1, doi: 10.1109/VETECF.2003.1285069.https://ieeexplore.ieee.org/document/1285069 

该部分将基于上述六种算法的运算结果来进行模拟仿真,由此来对 MIMO 系统中的各个算法性能进行比较分析。此处主要使用 MATLAB 对各个算法进行可视化仿真,主要的指标为不同信噪比情况下的误码率,运行速度,在不同传输天线数量情况下比较误码率,运行速度,由此来估算其实际传输过程之中的质量。我们控制

This part will carry out the simulation based on the operation results of the above six algorithms, so as to compare and analyze the performance of each algorithm in the MIMO system. Here we mainly use MATLAB to carry out visual simulation of each algorithm. The main indicators are the bit error rate and running speed under different signal-to-noise ratios, and compare the bit error rate and running speed under different numbers of transmission antennas. quality during actual transmission. we control

• 第一次仿真为 nT=4,nR=6, 传输的字长为 10^6bit 的情况下不同算法不同性能

• The first simulation is nT=4, nR=6, the transmitted word length is 10^6bit, the performance of different algorithms is different

• 第二次为 nT=6,nR=8, 传输的字长为 10^6bit 的情况下不同算法不同性能

• When the word length of the transmission is 10^6bit for the third time, the same algorithm remains unchanged under nT=4, and the bit error rate under nR=6-10 is different.

• 第三次为传输的字长为 10^6bit 的情况下,相同算法在 nT=4 不变,nR=6-10 之下的误码率的不同,

• The second time is nT=6, nR=8, and different algorithms have different performances when the transmitted word length is 106bit
在进行仿真模拟时,使用的信道为符合瑞利衰落的信道,传输的字长为 10^6bit,编码方式为QPSK 编码。

In the simulation, the channel used is the channel conforming to Rayleigh fading, the transmitted word length is 106bit, and the encoding method is QPSK encoding.

• First simulation
• 第一次仿真

![image](https://user-images.githubusercontent.com/95970863/199699190-bfae2a3d-d9a9-4721-9911-ce604cc90412.png)
![image](https://user-images.githubusercontent.com/95970863/199699228-a32a0272-56ee-4fd7-9a0a-ead2ee0d3f95.png)

以上两图,在 nT=4,nR=6, 传输的字长为 10^6bit 的情况下,VBLAST 算法的误 码率降低速度是最快的,SQRD 其次且和 VBLAST 相差不大,QRD 稍逊于 SQRD 但也差 不了太多。而使用 MMSE 检测法和 ZF 检测的线性算法总体上讲性能较差,包括 MMSEQR 在内,误码率明显高于前三种。

In the above two figures, in the case of nT=4, nR=6, and the transmitted word length is 10^6bit, VBLAST algorithm reduces the bit error rate the fastest, SQRD is second and is not much different from VBLAST, and QRD is slightly inferior Same as SQRD but not too much. However, the linear algorithm using MMSE detection method and ZF detection method has poor performance on the whole, including MMSEQR, the bit error rate is significantly higher than the first three.

![image](https://user-images.githubusercontent.com/95970863/199699445-d2718829-d405-4c5a-a8e8-4fc7d5fcfcd9.png)

上图显示,在 nT=6,nR=8, 传输的字长为 10^6bit 的情况下,VBLAST 算法的误码率降低速度 是最快的,SQRD 其次且和 VBLAST 相差不大,QRD 稍逊于 SQRD 但也差不了太多. 而天 线数量的增加使得误码率下降速度也更快

In the above figure it shows that in the case of nT=6, nR=8, and the transmitted word length is 10^6 bits, VBLAST algorithm reduces the bit error rate the fastest, SQRD is second and is not much different from VBLAST, and QRD is slightly inferior to SQRD But it's not too bad. The increase in the number of antennas makes the bit error rate drop faster

![image](https://user-images.githubusercontent.com/95970863/199699854-ff8ecfdb-d9dc-43ba-b4f7-f8be72e355ea.png)

上图展示了传输的字长为 10^6bit 的情况下,相同算法在 nT=4 不变,nR=6-10 之下的误码率 误码率与接受天线 nR 的数量的关系。当接收天线 nR 数量增加时,得到的信号的误码率更 低。这说明接收天线多有助于降低误码率

The figure above shows the relationship between the bit error rate and the number of receiving antennas nR when the word length of the transmission is 10^6bit, the same algorithm is unchanged at nT=4, and the bit error rate under nR=6-10. When the number of receiving antenna nR increases, the bit error rate of the obtained signal is lower. This shows that more receiving antennas help to reduce the bit error rate

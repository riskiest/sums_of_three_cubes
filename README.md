# sums_of_three_cubes

A non-official python implement of paper ["CRACKING THE PROBLEM WITH 33"](https://people.maths.bris.ac.uk/~maarb/papers/cubesv1.pdf) which present a computational solution to sums of three cubes problems for 33.

Now this project is just a toy python implement to the paper above. However, it accomplishes the gist. Improvement may or may not plan to be done.

现在的实现还远不够完美，基本还处于玩具阶段，但基本实现了大体解题思路，还有一些没有实现的部分如下：

1. [Montgomery‘s trick](https://en.wikipedia.org/wiki/Montgomery_modular_multiplication),这个我怀疑python的pow函数就有内置这个算法;
2. ~~M'; 由于M'的选用有很大的任意性，论文并没有给出M'的选用方法，所以这部分暂时没有去实现;~~M'只能选较小的素数，这部分主要在于判别式的计算资源与滤过的数值之间的平衡；
3. ~~算法根据d与$\sqrt {\alpha B}$的关系分成2部分;这样的话大于$\sqrt {\alpha B}$的素数的三次方表格就不需要计算了;~~这部分已经实现
4. C实现;这个不打算做，主要是python写的惯;
5. 多进程;程序结构本质是多进程的，有空弄下;
6. 日志;现在居然一点日志都没有
7. ~~表格的组成;现在很多都弄成表格，速度快了，但内存占用太大，计算不了太大的值;因此可能考虑放开一些表格，实时但重复计算，做时间-空间交易;~~现在已经能计算33的三次方根了

大家可以用http://www.asahi-net.or.jp/~KC2H-MSM/mathland/math04/matb0100.htm 找一些测试案例，但请注意，论文的方法只对$k\equiv3(mod\;9)$成立;

参考文献：

求解z^3=k(mod p)https://eprint.iacr.org/2009/457.pdf 以及其引用[3]

Rolandb在https://stackoverflow.com/questions/6752374/cube-root-modulo-p-how-do-i-do-this 的分享

一些代码引用修改自https://rosettacode.org/

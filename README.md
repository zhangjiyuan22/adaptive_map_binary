# Adaptive Map Binary Grid Search 使用说明
<br>
<br>
<br>
<br>
  总共有三个文件，包括生成map时使用的adaptive_map_generator.py，以及利用map做grid search时的chi2_calculator.cpp和binary_grid_ob16****.py，
另外还需要一个compile_all来编译chi2_calculator.cpp；此外还有一些测试用代码；
<br>
<br>
  另外还在项目中放了一个非常小的 map set（180张）以及一个事件的数据（ob161195），其grid search结果也已经放到result目录下，供测试用；
<br>
<br>
<br>

## 生成map：


#### 生成map时需要先mkdir map_set_**，然后在adaptive_map_generator.py中修改以下几个地方：


第15行处的box_size = ** ，每张map的大小是(2*box_size)*(2*box_size)，map的中心默认为binary的质心（在s和q均较大的情况下，也可以通过246/247行进行修改，将map中心调整到放大中心），目前同一次生成的一系列map还都是使用同样的box_size（后续会改进使得每张map的box_size自动选择），所以需要根据所选的log s、log q形成的caustic大小
来确保map的大小能装下所有magnification anomaly的区域；一般可以选用0.4/1.0/3.5等值；map之外用假设位于map中心的single近似放大倍数（s和q均较大,e.g.s=10,q=1,时需要考虑到map的长度单位是binary总质量对应的Einstein半径，而此时single计算公式的长度单位应是binary中一个质量对应的Einstein半径，所以需要进行转换）；此外也可以选择生成较大的map使得t_begin和t_end之间的数据全部在map里；
<br>

第234行处的if current_layer == ** : ，这里输入的数字代表map被加密到第**层之后强行停止加密；判断是否加密的标准是放大倍数的插值误差小于
<p align="center">
0.01*sqrt(μ) ,
</p> 
其中μ为当前位置的放大倍数，所以为了防止在靠近放大中心的位置不断加密，需要加密一定层之后强行停止的机制；加密到第11层就代表最小的小格边长为 
<p align="center">
2*box_size/(2^11),
</p> 
如果box_size = 0.4 ，那么current_layer == 11所实现的最小小格边长约为4*10^(-4)，这对于high-magnification event一般也是够用的，对普通的event可以自行选取这一值；
<br>
<br>
  第21行处的map_set_name = ‘'，需要修改为预先mkdir的那个准备存map的目录；
<br>
<br>
  第288-328行saveparm函数中修改log s / log q / log rho的范围和步长，也可以限制grid区域为特殊形状；
<br>
<br>
  第346行设置processes数目；
<br>

#### 然后直接python3 adaptive_map_generator.py即可运行，注意要这个python程序需要import MulensModel包

<br>
<br>

## 利用map做grid search：


#### chi2_calculator.cpp在正常使用中无需修改，也无需重复编译，只要在首次使用前通过./compile_all编译一次即可；注：compile_all中的路径需要修改；

#### 需要在binary_grid_ob16****.py中修改：
  
  第12行处设置所使用的map的目录名称；
  
  第126-129行处设置所使用的alpha初值；
  
  第176-186行以及第191/192行设置所使用的数据；
  
  第188/189/190行设置所使用的t0/u0/tE初值；
  
  第213-218行设置log s / log q / log rho的范围和步长（后面代码取其与所使用的map目录的交集）；
  
  第319行设置processes数目；
  
  第323行设置grid search结果存储位置及名称；

#### 然后直接python3 binary_grid_ob16****.py即可运行

#### 然后可以python3 draw_grid_search_result.py画图分析grid search的结果（需修改其12/13/15/16/137-140行）

<br>
<br>
<br>
<br>
<br>
<br>


# 另外，还有一些测试用的代码：

<br>
<br>
<br>

## 利用map插值得到放大倍数：

#### 在get_magnification.py中修改：
  
  第18/19行处设置使用哪张map进行插值；

  第13/14行处设置希望得到哪个(x,y)点的放大倍数；

#### 然后直接python3 get_magnification.py即可运行

<br>
<br>
<br>

## 将生成好的map可视化：

#### 在visualizing_adaptive_map.py 中修改：
  
  第45/46/48行处设置要可视化哪张map ；

#### 然后直接python3 visualizing_adaptive_map.py即可运行






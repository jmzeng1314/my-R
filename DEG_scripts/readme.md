## 一步法差异分析

有很多种选择，主要是继承自上面的normalization方法，一般来说挑选好了normalization方法就决定了选取何种差异分析方法，也并不强求弄懂统计学原理，它们都被包装到了对应的R包里面，主要是对R包的学习。

- edgeR (Robinson et al., 2010)
- DESeq / qDESeq2 (Anders and Huber, 2010, 2014)
- DEXSeq (Anders et al., 2012)
- limmaVoom
- Cuffdiff / Cuffdiff2 (Trapnell et al., 2013)
- PoissonSeq
- baySeq

## 制作表达矩阵

|                 | SRR1039508 | SRR1039509 | SRR1039512 | SRR1039513 | SRR1039516 | SRR1039517 | SRR1039520 | SRR1039521 |
| --------------- | ---------- | ---------- | ---------- | ---------- | ---------- | ---------- | ---------- | ---------- |
| ENSG00000000003 | 679        | 448        | 873        | 408        | 1138       | 1047       | 770        | 572        |
| ENSG00000000005 | 0          | 0          | 0          | 0          | 0          | 0          | 0          | 0          |
| ENSG00000000419 | 467        | 515        | 621        | 365        | 587        | 799        | 417        | 508        |
| ENSG00000000457 | 260        | 211        | 263        | 164        | 245        | 331        | 233        | 229        |
| ENSG00000000460 | 60         | 55         | 40         | 35         | 78         | 63         | 76         | 60         |
| ENSG00000000938 | 0          | 0          | 2          | 0          | 1          | 0          | 0          | 0          |
| ENSG00000000971 | 3251       | 3679       | 6177       | 4252       | 6721       | 11027      | 5176       | 7995       |
| ENSG00000001036 | 1433       | 1062       | 1733       | 881        | 1424       | 1439       | 1359       | 1109       |
| ENSG00000001084 | 519        | 380        | 595        | 493        | 820        | 714        | 696        | 704        |
| ENSG00000001167 | 394        | 236        | 464        | 175        | 658        | 584        | 360        | 269        |
| ENSG00000001460 | 172        | 168        | 264        | 118        | 241        | 210        | 155        | 177        |
| ENSG00000001461 | 2112       | 1867       | 5137       | 2657       | 2735       | 2751       | 2467       | 2905       |
| ENSG00000001497 | 524        | 488        | 638        | 357        | 676        | 806        | 493        | 475        |
| ENSG00000001561 | 71         | 51         | 211        | 156        | 23         | 38         | 134        | 172        |



第一列是基因ID，后面的列是各个样本。其中第一行尤为注意，最开头是一个空格(了解R里面read.table函数原理)

## 制作分组矩阵

|            | dex   | SampleName | cell    |
| ---------- | ----- | ---------- | ------- |
| SRR1039508 | untrt | GSM1275862 | N61311  |
| SRR1039509 | trt   | GSM1275863 | N61311  |
| SRR1039512 | untrt | GSM1275866 | N052611 |
| SRR1039513 | trt   | GSM1275867 | N052611 |
| SRR1039516 | untrt | GSM1275870 | N080611 |
| SRR1039517 | trt   | GSM1275871 | N080611 |
| SRR1039520 | untrt | GSM1275874 | N061011 |
| SRR1039521 | trt   | GSM1275875 | N061011 |

记住要跟上面的表达矩阵的样本名对应！！！

只有第一列是需要看的，其余的无所谓。

根据分组信息，是需要自己指定比对信息的，比如上面的分组矩阵，需要指定 `-c 'trt-untrt'`

## 下载差异分析脚本

```shell
wget  https://github.com/jmzeng1314/my-R/blob/master/DEG_scripts/run_DEG.R
wget  https://github.com/jmzeng1314/my-R/blob/master/DEG_scripts/tair/exprSet.txt
wget  https://github.com/jmzeng1314/my-R/blob/master/DEG_scripts/tair/group_info.txt
Rscript ../run_DEG.R -e exprSet.txt -g group_info.txt -c 'Day1-Day0' -s counts  -m DESeq2
```

如果是转录组的raw counts数据，就选择 -s counts，如果是芯片等normalization好的表达矩阵数据，用默认参数即可。下面是例子：

```shell

# Rscript run_DEG.R -e airway.expression.txt -g airway.group.txt -c 'trt-untrt' -s counts -m DESeq2
# Rscript run_DEG.R -e airway.expression.txt -g airway.group.txt -c 'trt-untrt' -s counts -m edgeR
# Rscript run_DEG.R -e sCLLex.expression.txt -g sCLLex.group.txt -c 'progres.-stable'
# Rscript run_DEG.R -e sCLLex.expression.txt -g sCLLex.group.txt -c 'progres.-stable' -m t.test
```

对于转录组的raw counts数据，有DEseq2包和edgeR包可供选择。对于芯片等normalization好的表达矩阵数据，有limma和t.test供选择。



## 重要的脚本

比如 `create_testData.R` 里面有如何得到表达矩阵和分组矩阵的内容。


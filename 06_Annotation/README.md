# 注释流程使用方法

## toy data 测试

准备了马铃薯的单条染色体和一个转录组数据进行测试，由于braker本身的限制无法在测试数据中运行，所以默认提供的流程不进行braker，
只运行maker的第一轮测试

```bash
cp -r /public/agis/huangsanwen_group/baozhigui/Pipeline/Annotation ./
cd Annotation

## 修改配置文件
sed -i "s|^outdir:|outdir: \"`pwd`\"|g" config.yaml

## 加载snakemake环境
source /public/home/baozhigui/miniconda3/bin/activate snakemake

## dry run测试任务是否逻辑正常，出现任务数统计即可正常运行
snakemake -npr -s Snakefile

## 运行-j 指定最大线程
nohup snakemake -s Snakefile -j {threads} --stats snakejob.stats >&2 2>> snakejob.log

```
 

## 实际项目测试

### 拷贝流程

```bash
cp -r /public/agis/huangsanwen_group/baozhigui/Pipeline/Annotation ./
cd Annotation

## 修改配置文件
sed -i "s|^outdir:|outdir: \"`pwd`\"|g" config.yaml
```




### 准备数据

（以下均为相对路径,以DM为例，链接到相应目录，没有的话请新建目录）
- 基因组文件: `DM.fa`(必须fa后缀)
- 二代转录组原始数据: `rawdata/04.rna-seq/tissue1.R1.fastq.gz`,`rawdata/04.rna-seq/tissue1.R2.fastq.gz`(支持多个组织)
- 同源物种蛋白: `rawdata/05.Annotation/proteins.fasta`(所有蛋白存放为一个文件即可)

### 填写配置文件 config.yaml

```yaml
rawdata:
    ## rawdata/04.rna-seq/{sample}.R1.fastq.gz
    ## 二代转录组
    RNA:
        - "potato_stolon1"
        - "potato_stolon2"

# 基因组序列前缀和输出前缀
name: "DM"

# 参数
parameters:
    EDTA:
        ## final step [all,filter,final,anno]
        step: "all"
        ## rice,maize,others
        # 物种，除水稻玉米以外都是others
        species: "others"
        # 0表示不跑repeatmodeler，1表示运行
        sensitive: 0

        ## 准确的已知TE库
        lib:
        ## 已注释的cds序列
        cds:
        # lib: "maizeTE02052020"
        # cds: "Zea_mays.B73_RefGen_v4.cdna.all.fa"
```
### 修改snakefile

去除以下两行的注释，就会运行braker和maker的round2
```bash
        #"results/03.annotation/05.braker/augustus.hints.gff3",
        expand("results/03.annotation/06.maker/round1/{name}.maker.output/{name}.maker.gff",name=name),
        #expand("results/03.annotation/06.maker/round2/{name}.maker.output/{name}.maker.gff",name=name)
```

### 运行

```bash
## 加载snakemake环境
source /public/home/baozhigui/miniconda3/bin/activate snakemake

## dry run测试任务是否逻辑正常，出现任务数统计即可正常运行
snakemake -npr -s Snakefile

## 运行-j 指定最大线程
nohup snakemake -s Snakefile -j {threads} --stats snakejob.stats >&2 2>> snakejob.log
```





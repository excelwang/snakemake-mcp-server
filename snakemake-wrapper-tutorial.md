# 教程：将你的分析工具制作成一个可复用的 Snakemake 包装器

## 第一部分：为什么要做这件事？

你可能有一个非常强大的脚本（例如，一个用 Python、R 或其他语言写的脚本），它能出色地完成一项特定任务。但现在，你希望：

1.  **分享给同事**：让他们能轻松地在自己的数据上使用你的工具，而不需要修改你的脚本代码。
2.  **重复使用**：在不同的项目和不同的数据上，快速地调用这个工具。
3.  **自动化**：将这个工具作为大型分析流程中的一个步骤，实现从原始数据到最终结果的全程自动化。

“包装器”（Wrapper）就是为了实现这些目标而设计的“标准盒子”。它将你的脚本、依赖环境和说明书打包在一起，让一个叫做 **Snakemake** 的工作流管理器能够理解并运行它。

你不需要成为 Snakemake 专家，只需要遵循以下步骤，就能为你的工具打造一个这样的“标准盒子”。

---

## 第二部分：准备工作——一个简单的“工具脚本”

假设你有一个 Python 脚本 `my_tool.py`，它的功能是：读取一个 CSV 文件，筛选出某一列大于特定阈值的数据，然后保存为新的 CSV 文件。

这个脚本目前可能长这样（**“封装前”的样子**）：

```python
# my_tool.py - 封装前的样子
import pandas as pd

# 问题：输入/输出文件和参数都是“写死的”
input_file = "/path/to/my/data.csv"
output_file = "/path/to/my/result.csv"
threshold = 50

# 核心逻辑
df = pd.read_csv(input_file)
filtered_df = df[df["value_column"] > threshold]
filtered_df.to_csv(output_file, index=False)

print(f"处理完成！结果保存在 {output_file}")
```

这个脚本很难直接给别人用，因为他们必须手动修改脚本里的文件路径和 `threshold` 值。我们的目标就是改造它。

---

## 第三部分：包装器的标准结构

一个标准的包装器是一个目录，它包含以下部分（就像我们看到的 `rasterio/clip` 目录一样）：

```
my_awesome_tool/
├── environment.yaml   # 软件“购物清单”
├── meta.yaml          # “说明书”
├── wrapper.py         # 改造后的“工具脚本”
└── test/              # “测试车间”
    ├── Snakefile      # 测试流程脚本
    └── data.csv       # 测试用的输入数据
```

**现在，我们开始动手创建！**

---

## 第四部分：一步步创建你的包装器

### 第 1 步：创建目录和文件

首先，创建上述的目录结构和空的 `environment.yaml`, `meta.yaml`, `wrapper.py` 文件。

### 第 2 步：定义环境 (`environment.yaml`)

这个文件告诉 Snakemake 你的脚本需要哪些软件才能运行。这保证了任何人在任何电脑上都能配置出完全相同的运行环境。

对于我们这个简单的例子，脚本用到了 `pandas` 库。所以 `environment.yaml` 文件内容如下：

```yaml
# environment.yaml
channels:
  - conda-forge
dependencies:
  - pandas =2.2.2 # 建议指定版本以保证稳定性
```

**类比 `rasterio/clip`**：它的 `environment.yaml` 文件里列出了 `rasterio`, `geopandas` 等库，因为它的 `wrapper.py` 脚本需要这些库来工作。

### 第 3 步：改造核心脚本 (`wrapper.py`)

这是最关键的一步。我们将移除所有“写死”的路径和参数，让脚本从 Snakemake 那里接收它们。

Snakemake 会在运行时创建一个特殊的全局对象 `snakemake`，你的脚本可以通过它来访问一切外部信息。

**“封装后”的 `wrapper.py` 应该长这样：**

```python
# wrapper.py - 封装后的样子
import pandas as pd

# 1. 从 snakemake 对象获取输入文件路径
# snakemake.input 是一个列表，[0] 代表第一个输入文件
input_file = snakemake.input[0] 

# 2. 从 snakemake 对象获取输出文件路径
# snakemake.output 是一个列表，[0] 代表第一个输出文件
output_file = snakemake.output[0]

# 3. 从 snakemake 对象获取参数
# snakemake.params 是一个字典，可以通过名字访问
threshold = snakemake.params.threshold

# 4. (可选，但推荐) 将所有打印信息重定向到日志文件
# 这能保持工作流输出的整洁
log_file = snakemake.log[0]
with open(log_file, "w") as f:
    # 将脚本的输出重定向到日志文件
    import sys
    sys.stderr = sys.stdout = f

    # --- 核心逻辑 (和原来一样，无需改动) ---
    df = pd.read_csv(input_file)
    filtered_df = df[df["value_column"] > threshold]
    filtered_df.to_csv(output_file, index=False)

    print(f"处理完成！输入: {input_file}, 输出: {output_file}, 阈值: {threshold}")
```

**与外部世界交互的桥梁——`snakemake` 对象详解：**

*   `snakemake.input`: 一个文件路径列表或字典。它由 Snakemake 根据工作流的规则自动提供。
    *   `snakemake.input[0]`: 获取第一个输入文件。
    *   `snakemake.input.my_data`: 如果输入在规则中被命名了（如 `my_data="path/to/file"`），可以用名字访问。
*   `snakemake.output`: 与 `snakemake.input` 类似，用于获取输出文件的路径。你的脚本**必须**生成这个路径对应的文件。
*   `snakemake.params`: 一个字典，用于接收**非文件类**的参数，如数字、字符串等。
    *   `snakemake.params.threshold`: 获取名为 `threshold` 的参数。
*   `snakemake.log`: 日志文件的路径。将脚本的屏幕输出（print、错误信息等）写入这里是最佳实践。
*   `snakemake.threads`: 如果你的工具能利用多线程，可以通过它获取 Snakemake 分配的 CPU 核心数。

**总结：你的脚本不再关心文件具体在哪里，或者参数值到底是多少。它只跟 `snakemake` 这个“经纪人”对话，由“经纪人”负责传递所有信息。**

### 第 4 步：撰写说明书 (`meta.yaml`)

这个文件用统一的格式描述你的包装器是做什么的，需要什么，产出什么。这对于使用者来说至关重要。

```yaml
# meta.yaml
name: "My Awesome Tool"
description: "一个筛选 CSV 文件中数值的工具。"
authors:
  - "你的名字"

# 描述需要哪些输入文件
input:
  - "一个包含'value_column'列的 CSV 文件。"

# 描述会生成哪些输出文件
output:
  - "一个只包含筛选后数据的 CSV 文件。"

# 描述可以接受哪些额外参数
params:
  - threshold: "用于筛选'value_column'的数值阈值 (默认: 50)。"
```

**类比 `rasterio/clip`**：它的 `meta.yaml` 清晰地描述了输入（栅格、矢量等）、输出（裁剪后的TIFF）和各种参数（边界、缓冲等）。

### 第 5 步：创建测试 (`test/` 目录)

测试是保证你的包装器能正常工作的关键。

1.  **准备测试数据**：在 `test/` 目录下，创建一个简单的 `data.csv` 文件。
    ```csv
    id,value_column
    A,10
    B,80
    C,45
    D,120
    ```

2.  **编写测试 `Snakefile`**：`Snakefile` 是 Snakemake 的“流程定义文件”。在 `test/` 目录下创建它，内容如下：

    ```python
    # test/Snakefile

    # 定义一个规则，命名为 test_my_awesome_tool
    rule test_my_awesome_tool:
        # 定义输入文件
        input:
            "data.csv"
        # 定义期望的输出文件
        output:
            "result.csv"
        # 定义要传递的参数
        params:
            threshold=50
        # 定义日志文件
        log:
            "test.log"
        # 关键！告诉 Snakemake 使用哪个包装器
        wrapper:
            "file:../" # "file:" 表示使用本地文件，"../" 指向包装器根目录
    ```

3.  **运行测试**：
    打开终端，进入 `test/` 目录，然后运行以下命令：

    ```bash
    # 这个命令会：
    # 1. --use-conda: 自动使用 conda 创建在 environment.yaml 中定义的环境
    # 2. --conda-frontend mamba: 使用更快的 mamba 来安装环境
    # 3. -c1: 使用一个核心来执行
    snakemake --use-conda --conda-frontend mamba -c1 
    ```

    如果一切顺利，Snakemake 会：
    1.  读取 `Snakefile`，看到 `test_my_awesome_tool` 规则。
    2.  发现需要 `environment.yaml` 中定义的环境，并自动创建它。
    3.  调用 `wrapper.py` 脚本，并将 `input`, `output`, `params` 等信息打包到 `snakemake` 对象中传递给它。
    4.  脚本成功执行后，你会在 `test/` 目录下看到 `result.csv` 和 `test.log` 文件。

    `result.csv` 的内容应该是：
    ```csv
    id,value_column
    B,80
    D,120
    ```

---

## 总结

恭喜你！你已经成功地将一个独立的脚本封装成了一个结构标准、易于分享和可被自动流程调用的 Snakemake 包装器。

**核心要点回顾：**

1.  **标准化结构**：使用 `wrapper.py`, `environment.yaml`, `meta.yaml`, `test/` 的标准目录结构。
2.  **解耦**：在 `wrapper.py` 中，通过 `snakemake` 对象（`snakemake.input`, `snakemake.output`, `snakemake.params`）来与外部世界通信，彻底告别写死的路径和参数。
3.  **环境隔离**：在 `environment.yaml` 中明确声明所有软件依赖，保证可复现性。
4.  **文档清晰**：在 `meta.yaml` 中为你的用户提供清晰的“说明书”。
5.  **测试驱动**：提供一个 `test/` 目录，包含测试数据和 `Snakefile`，确保你的包装器稳定可靠。

---

## 附录：Snakemake 包装器中的 `snakemake` 对象详解

这份附录是 `snakemake` 对象的速查手册，它是在 `wrapper.py` 脚本中与工作流交互的唯一桥梁。

### `snakemake.input`

*   **类型**: `snakemake.io.InputFiles` (行为类似列表和字典的混合体)
*   **描述**: 包含所有输入文件的路径。
*   **用法与示例**:
    *   **按索引访问** (当输入是匿名的):
        *   `Snakefile`:
            ```python
            rule my_rule:
                input: "data/a.txt", "data/b.txt"
                output: "results/c.txt"
                wrapper: "file:path/to/wrapper"
            ```
        *   `wrapper.py`:
            ```python
            first_input = snakemake.input[0]  # "data/a.txt"
            second_input = snakemake.input[1] # "data/b.txt"
            # 也可以迭代
            for f in snakemake.input:
                print(f)
            ```
    *   **按名称访问** (当输入被命名):
        *   `Snakefile`:
            ```python
            rule my_rule:
                input:
                    reads="data/reads.fastq",
                    reference="data/genome.fasta"
                output: "results/aln.bam"
                wrapper: "file:path/to/wrapper"
            ```
        *   `wrapper.py`:
            ```python
            reads_file = snakemake.input.reads      # "data/reads.fastq"
            ref_file = snakemake.input.reference  # "data/genome.fasta"
            ```

### `snakemake.output`

*   **类型**: `snakemake.io.OutputFiles` (行为类似列表和字典的混合体)
*   **描述**: 包含所有输出文件的路径。你的脚本**必须**创建这些文件。
*   **用法与示例**: 用法与 `snakemake.input` 完全相同。
    *   **按索引访问**:
        *   `Snakefile`: `output: "results/plot.png", "results/table.csv"`
        *   `wrapper.py`: `plot_path = snakemake.output[0]`
    *   **按名称访问**:
        *   `Snakefile`: `output: plot="results/plot.png", table="results/table.csv"`
        *   `wrapper.py`: `plot_path = snakemake.output.plot`

### `snakemake.params`

*   **类型**: `snakemake.io.Params` (行为类似字典)
*   **描述**: 包含所有在 `params` 块中定义的**非文件类**参数。
*   **用法与示例**:
    *   `Snakefile`:
        ```python
        rule my_rule:
            input: "data.csv"
            output: "result.csv"
            params:
                threshold=50,
                method="linear",
                columns=["col_a", "col_b"]
            wrapper: "file:path/to/wrapper"
        ```
    *   `wrapper.py`:
        ```python
        threshold = snakemake.params.threshold  # 50
        method = snakemake.params.method      # "linear"
        cols = snakemake.params.columns       # ["col_a", "col_b"]
        ```

### `snakemake.wildcards`

*   **类型**: `snakemake.io.Wildcards` (行为类似字典)
*   **描述**: 包含从输入/输出文件路径中匹配到的**通配符**的值。这是实现流程自动扩展的关键。
*   **用法与示例**:
    *   `Snakefile`:
        ```python
        # {sample} 和 {ext} 就是通配符
        rule my_rule:
            input: "data/{sample}.{ext}"
            output: "results/{sample}/report.html"
            wrapper: "file:path/to/wrapper"
        ```
        当 Snakemake 需要创建 `results/sample1/report.html` 时：
    *   `wrapper.py`:
        ```python
        sample_name = snakemake.wildcards.sample  # "sample1"
        extension = snakemake.wildcards.ext       # "txt" (假设输入是 data/sample1.txt)
        print(f"Processing sample: {sample_name}")
        ```

### `snakemake.log`

*   **类型**: `snakemake.io.Log` (行为类似列表和字典的混合体)
*   **描述**: 包含日志文件的路径。强烈建议将所有屏幕输出（标准输出/错误）重定向到此文件。
*   **用法与示例**:
    *   `Snakefile`: `log: "logs/my_rule.log"`
    *   `wrapper.py`:
        ```python
        import sys
        log_file = snakemake.log[0]
        with open(log_file, "w") as log:
            sys.stderr = sys.stdout = log
            # 现在所有的 print 和错误信息都会写入 logs/my_rule.log
            print("This goes to the log file.")
        ```

### `snakemake.threads`

*   **类型**: `int`
*   **描述**: 当前任务被分配到的 CPU 线程数。
*   **用法与示例**:
    *   `Snakefile`:
        ```python
        rule align_reads:
            ...
            threads: 8
            wrapper: "file:path/to/wrapper"
        ```
    *   `wrapper.py`:
        ```python
        from snakemake.shell import shell
        # 将线程数传递给命令行工具
        shell(f"bwa_tool -t {snakemake.threads} ...") 
        ```

### `snakemake.resources`

*   **类型**: `snakemake.io.Resources` (行为类似字典)
*   **描述**: 访问为任务定义的任意资源，如内存、GPU等。
*   **用法与示例**:
    *   `Snakefile`:
        ```python
        rule deep_learning:
            ...
            resources:
                mem_mb=8000,
                gpu=1
            wrapper: "file:path/to/wrapper"
        ```
    *   `wrapper.py`:
        ```python
        mem = snakemake.resources.mem_mb  # 8000
        gpus = snakemake.resources.gpu    # 1
        print(f"Requesting {mem}MB of memory.")
        ```

### `snakemake.config`

*   **类型**: `dict`
*   **描述**: 访问整个工作流的全局配置字典。通常通过 `snakemake --configfile config.yaml` 载入。
*   **用法与示例**:
    *   `config.yaml`:
        ```yaml
        global_threshold: 0.95
        samples_dir: "data/samples"
        ```
    *   `Snakefile`: (无需特殊声明，自动可用)
    *   `wrapper.py`:
        ```python
        # 访问配置
        threshold = snakemake.config["global_threshold"] # 0.95
        samples_path = snakemake.config["samples_dir"]   # "data/samples"
        ```

### `snakemake.rule`

*   **类型**: `str`
*   **描述**: 正在执行的规则的名称。
*   **用法与示例**:
    *   `Snakefile`: `rule my_awesome_rule:`
    *   `wrapper.py`: `print(f"Executing rule: {snakemake.rule}")` # "my_awesome_rule"

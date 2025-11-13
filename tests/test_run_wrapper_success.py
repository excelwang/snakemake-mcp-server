import pytest
import os
import shutil

@pytest.mark.asyncio
async def test_run_wrapper_success(test_files, run_wrapper_test):
    """测试通过直接函数调用成功执行wrapper"""
    wrapper_runner, workdir = run_wrapper_test
    
    # 复制测试输入文件到 run_wrapper_test 创建的临时目录
    input_filename = os.path.basename(test_files['input'])
    output_filename = os.path.basename(test_files['output'])
    input_path = os.path.join(workdir, input_filename)
    output_path = os.path.join(workdir, output_filename)
    
    shutil.copy2(test_files['input'], input_path)

    result = await wrapper_runner(
        wrapper_name="bio/samtools/faidx",
        inputs=[input_filename],
        outputs=[output_path],
        params={}
    )

    # 验证结果
    assert 'status' in result, "Result should have status attribute"

    # 驗證執行狀態
    assert result['status'] == 'success', \
        f"Expected success, got {result.get('status')}: {result.get('error_message')}"

    # 驗證輸出文件
    assert os.path.exists(output_path), \
        f"Output file should be created: {output_path}"

    # 驗證文件內容
    with open(output_path, 'r') as f:
        content = f.read().strip()
        assert len(content) > 0, "Output file should not be empty"
        assert '\t' in content, "FAI file should be tab-delimited"

import pytest
import os
import shutil

@pytest.mark.asyncio
async def test_run_wrapper_error_handling(test_files, run_wrapper_test):
    """测试直接函数调用的错误处理"""
    wrapper_runner, workdir = run_wrapper_test
    
    # 复制测试输入文件到 run_wrapper_test 创建的临时目录
    input_filename = os.path.basename(test_files['input'])
    output_filename = os.path.basename(test_files['output'])
    input_path = os.path.join(workdir, input_filename)
    output_path = os.path.join(workdir, output_filename)
    
    shutil.copy2(test_files['input'], input_path)

    # The function should return a failure result rather than raising an exception
    result = await wrapper_runner(
        wrapper_name="",  # 无效参数
        inputs=[input_filename],
        outputs=[output_path],
        params={}
    )

    # Should return failed status
    assert result["status"] == "failed"
    assert "wrapper_name must be a non-empty string" in str(result.get("error_message", ""))
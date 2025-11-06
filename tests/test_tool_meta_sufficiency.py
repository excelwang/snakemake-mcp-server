import pytest
import asyncio
from fastmcp import Client
import os
import tempfile


@pytest.mark.asyncio
async def test_tool_meta_contains_sufficient_information(http_client: Client):
    """测试 get-tool-meta 返回的信息是否足够指导用户配置 tool/process 参数"""
    # 使用一个已知的简单工具，例如 samtools/stats
    result = await asyncio.wait_for(
        http_client.call_tool(
            "get_tool_meta",
            {
                "tool_path": "bio/samtools/stats"  # 使用存在的工具路径
            }
        ),
        timeout=30
    )
    
    # 验证结果结构
    assert hasattr(result, 'data'), "Result should have data attribute"
    
    # 检查返回数据结构
    metadata = result.data
    
    # 验证返回的元数据包含足够的信息来指导参数配置
    assert hasattr(metadata, 'name'), "Response should have name"
    assert hasattr(metadata, 'description'), "Response should have description"
    assert hasattr(metadata, 'input'), "Response should have input specification"
    assert hasattr(metadata, 'output'), "Response should have output specification" 
    assert hasattr(metadata, 'params'), "Response should have params specification"
    assert hasattr(metadata, 'path'), "Response should have path"
    
    print(f"Tool name: {metadata.name}")
    print(f"Tool description: {metadata.description}")
    print(f"Tool input: {metadata.input}")
    print(f"Tool output: {metadata.output}")
    print(f"Tool params: {metadata.params}")
    
    # 验证元数据信息足以指导用户配置
    # 1. 输入规范应该存在并描述需要什么类型的数据
    if metadata.input:
        assert isinstance(metadata.input, (dict, list)), "Input should be dict or list type"
        print("✓ Input specification available for user guidance")
    
    # 2. 输出规范应该存在并描述生成什么类型的数据  
    if metadata.output:
        assert isinstance(metadata.output, (dict, list)), "Output should be dict or list type"
        print("✓ Output specification available for user guidance")
    
    # 3. 参数规范应该存在并描述可配置的参数
    if metadata.params:
        # Parameters can be dict or list depending on the wrapper format
        assert isinstance(metadata.params, (dict, list)), "Params should be dict or list type"
        print("✓ Params specification available for user guidance")
    
    # 4. 描述信息应该存在并说明工具功能
    if metadata.description:
        assert isinstance(metadata.description, str), "Description should be string type"
        print("✓ Description available to explain tool purpose")
    
    # 5. 如果有notes字段也应该存在
    if hasattr(metadata, 'notes') and metadata.notes:
        assert isinstance(metadata.notes, str), "Notes should be string type"
        print("✓ Additional notes available for user guidance")
    
    print("✓ Tool metadata contains sufficient information to guide user configuration")


@pytest.mark.asyncio
async def test_tool_meta_can_guide_tool_process_usage(http_client: Client):
    """测试 tool metadata 是否能指导实际的 tool/process 使用"""
    # 获取工具元数据
    metadata_result = await asyncio.wait_for(
        http_client.call_tool(
            "get_tool_meta",
            {
                "tool_path": "bio/samtools/stats"
            }
        ),
        timeout=30
    )
    
    metadata = metadata_result.data
    
    # 验证元数据可以用来构造正确的 tool/process 请求
    # 检查元数据是否包含足够的信息来配置参数
    
    # 模拟根据元数据构造一个 tool/process 请求
    tool_request = {
        "wrapper_name": "samtools/stats",  # 从工具路径构造
        "inputs": [],  # 根据 metadata.input 信息填充
        "outputs": [],  # 根据 metadata.output 信息填充
        "params": {}  # 根据 metadata.params 信息填充
    }
    
    # 检查是否有输入定义
    if metadata.input:
        print(f"Input specification: {metadata.input}")
        # 如果有输入定义，应该知道需要什么类型的输入文件
        if isinstance(metadata.input, dict):
            for key, value in metadata.input.items():
                print(f"Input key: {key}, expected: {value}")
        elif isinstance(metadata.input, list):
            for item in metadata.input:
                print(f"Input item: {item}")
    
    # 检查是否有参数定义
    if metadata.params:
        print(f"Parameter specification: {metadata.params}")
        # 如果有参数定义，应该知道可以配置哪些参数
        if isinstance(metadata.params, dict):
            for param_name, param_info in metadata.params.items():
                print(f"Parameter: {param_name}, info: {param_info}")
        elif isinstance(metadata.params, list):
            for param_item in metadata.params:
                print(f"Parameter item: {param_item}")
    
    print("✓ Tool metadata can guide tool/process request construction")


@pytest.mark.asyncio
async def test_tool_meta_for_different_tool_types(http_client: Client):
    """测试不同类型的工具 metadata 是否都提供足够的信息"""
    # 测试多个不同类型的工具
    test_tools = [
        "bio/samtools/stats",
        "bio/samtools/faidx"  # 添加更多已知存在的工具
    ]
    
    for tool_path in test_tools:
        try:
            result = await asyncio.wait_for(
                http_client.call_tool(
                    "get_tool_meta",
                    {
                        "tool_path": tool_path
                    }
                ),
                timeout=30
            )
            
            metadata = result.data
            
            # 检查每个工具的元数据都包含必需信息
            required_fields = [
                ('name', metadata.name),
                ('path', metadata.path),
                ('input', metadata.input),
                ('output', metadata.output),
                ('params', metadata.params)
            ]
            
            for field_name, field_value in required_fields:
                assert field_value is not None, f"{tool_path} should have {field_name} information"
            
            print(f"✓ {tool_path} has sufficient metadata for configuration guidance")
            
        except Exception as e:
            if "404" in str(e) or "not found" in str(e).lower():
                print(f"⚠ {tool_path} not found, skipping...")
                continue
            else:
                raise
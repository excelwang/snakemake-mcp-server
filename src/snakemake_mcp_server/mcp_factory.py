from .api.main import create_native_fastapi_app
from fastmcp import FastMCP

def create_mcp_from_fastapi(wrappers_path: str, workflows_dir: str):
    """
    Create an MCP server from the native FastAPI application.
    This follows the recommended pattern from FastMCP documentation.
    """
    # First create the native FastAPI app
    fastapi_app = create_native_fastapi_app(wrappers_path, workflows_dir)
    
    # Convert to MCP server
    mcp = FastMCP.from_fastapi(app=fastapi_app)
    
    return mcp

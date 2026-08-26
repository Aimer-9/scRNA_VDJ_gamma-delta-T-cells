#!/usr/bin/env bash
set -euo pipefail

# Keep the MCP endpoint URL and default chat_id outside this script because
# they contain secrets/project-specific IDs. The env file is ignored by git.
config_file="${WECOM_CONFIG:-$(dirname "${BASH_SOURCE[0]}")/env/wecom_mcp.env}"
if [[ -f "$config_file" ]]; then
    # shellcheck source=/dev/null
    source "$config_file"
fi

url="${WECOM_MCP_URL:-}"
mcp_tool="${WECOM_MCP_TOOL:-message_send}"
wecom_chat_id="${WECOM_CHAT_ID:-}"

# Some Streamable HTTP MCP servers return an Mcp-Session-Id header after
# initialize. Reuse it for subsequent tools/list and tools/call requests.
mcp_session_id="${WECOM_MCP_SESSION_ID:-}"

usage() {
    cat <<'USAGE'
Usage:
  bash scripts/wecom_mcp.sh send "message text"
  bash scripts/wecom_mcp.sh send --chat-id CHAT_ID "message text"
  bash scripts/wecom_mcp.sh list-tools
  bash scripts/wecom_mcp.sh list-groups [begin_time] [end_time]
  bash scripts/wecom_mcp.sh list-messages --chat-id CHAT_ID [begin_time] [end_time]
  bash scripts/wecom_mcp.sh "message text"

Environment:
  WECOM_CONFIG         Config file to source, default: scripts/env/wecom_mcp.env
  WECOM_MCP_URL       MCP Streamable HTTP URL
  WECOM_MCP_TOOL       MCP tool name used for sending, default: message_send
  WECOM_CHAT_ID        Required chat ID for receiving the message

Examples:
  bash scripts/wecom_mcp.sh list-tools
  bash scripts/wecom_mcp.sh list-groups
  bash scripts/wecom_mcp.sh list-messages --chat-id CHAT_ID_HERE
  WECOM_CHAT_ID='CHAT_ID_HERE' bash scripts/wecom_mcp.sh send "job finished"
USAGE
}

require_cmd() {
    command -v "$1" >/dev/null 2>&1 || {
        echo "[ERROR] Required command not found: $1" >&2
        return 1
    }
}

json_escape() {
    # Build JSON strings through Python so quotes, newlines, and non-ASCII
    # message text are escaped correctly.
    python3 -c 'import json, sys; print(json.dumps(sys.stdin.read(), ensure_ascii=False))'
}

extract_json_response() {
    # Streamable HTTP can return either normal JSON or server-sent events.
    # Normalize both forms to one JSON object so downstream grep/tests work.
    python3 -c '
import json
import sys

body = sys.stdin.read().strip()
if not body:
    raise SystemExit("empty MCP response")

if body.startswith("event:") or "\ndata:" in body or body.startswith("data:"):
    data = []
    for line in body.splitlines():
        if line.startswith("data:"):
            chunk = line[5:].strip()
            if chunk and chunk != "[DONE]":
                data.append(chunk)
    body = "\n".join(data).strip()

decoder = json.JSONDecoder()
obj, _ = decoder.raw_decode(body)
print(json.dumps(obj, ensure_ascii=False))
'
}

mcp_post() {
    # Send one JSON-RPC request to the MCP Streamable HTTP endpoint. Headers
    # must advertise JSON/SSE support; otherwise this endpoint returns
    # "Not Acceptable: Client must accept application/json".
    local payload="$1"
    local header_file body_file response_json

    [[ -n "$url" ]] || {
        echo "[ERROR] WECOM_MCP_URL is required. Put it in scripts/env/wecom_mcp.env or set WECOM_MCP_URL." >&2
        return 1
    }

    header_file="$(mktemp)"
    body_file="$(mktemp)"

    if [[ -n "$mcp_session_id" ]]; then
        curl -sS -D "$header_file" -o "$body_file" \
            "$url" \
            -H 'Accept: application/json, text/event-stream' \
            -H 'Content-Type: application/json' \
            -H "Mcp-Session-Id: ${mcp_session_id}" \
            --data-binary "$payload"
    else
        curl -sS -D "$header_file" -o "$body_file" \
            "$url" \
            -H 'Accept: application/json, text/event-stream' \
            -H 'Content-Type: application/json' \
            --data-binary "$payload"
    fi

    if [[ -z "$mcp_session_id" ]]; then
        # Header names are case-insensitive. awk extracts the session ID from
        # the first response that provides it.
        mcp_session_id="$(
            awk 'BEGIN{IGNORECASE=1} /^Mcp-Session-Id:/ {sub(/\r$/, "", $2); print $2; exit}' "$header_file"
        )"
    fi

    response_json="$(extract_json_response < "$body_file")"
    rm -f "$header_file" "$body_file"

    printf '%s\n' "$response_json"
    if printf '%s\n' "$response_json" | grep -q '"error"[[:space:]]*:'; then
        return 1
    fi
}

mcp_initialize() {
    # MCP requires initialize before tool calls. This script does not need any
    # special client capabilities for simple message/list operations.
    local payload

    payload='{
      "jsonrpc": "2.0",
      "id": 1,
      "method": "initialize",
      "params": {
        "protocolVersion": "2025-03-26",
        "capabilities": {},
        "clientInfo": {
          "name": "wecom_mcp.sh",
          "version": "1.0.0"
        }
      }
    }'
    mcp_post "$payload" >/dev/null
}

mcp_list_tools() {
    # Print the server's tools/list response. Use this when the send tool name
    # or schema changes upstream.
    local payload

    mcp_initialize
    payload='{"jsonrpc":"2.0","id":2,"method":"tools/list","params":{}}'
    mcp_post "$payload"
}

default_begin_time() {
    # WeCom chat history tools only support a recent time window. Default to
    # the last 24 hours so list-groups/list-messages work without arguments.
    date -d '1 day ago' '+%F %T'
}

default_end_time() {
    date '+%F %T'
}

mcp_call_tool() {
    # Generic MCP tools/call wrapper. arguments_json must already be a JSON
    # object string matching the selected tool's input schema.
    local tool="$1"
    local arguments_json="$2"
    local payload

    payload="$(
        TOOL="$tool" ARGUMENTS_JSON="$arguments_json" python3 -c '
import json
import os

payload = {
    "jsonrpc": "2.0",
    "id": 3,
    "method": "tools/call",
    "params": {
        "name": os.environ["TOOL"],
        "arguments": json.loads(os.environ["ARGUMENTS_JSON"])
    }
}
print(json.dumps(payload, ensure_ascii=False))
'
    )"

    mcp_initialize
    mcp_post "$payload"
}

mcp_list_groups() {
    # Call chat_groups_list to find recent group conversation IDs. The visible
    # group name is not a valid chat_id; use the returned conversation ID.
    local begin_time="${1:-$(default_begin_time)}"
    local end_time="${2:-$(default_end_time)}"
    local args_json

    require_cmd curl
    require_cmd python3

    args_json="$(
        BEGIN_TIME="$begin_time" END_TIME="$end_time" python3 -c '
import json
import os
print(json.dumps({
    "begin_time": os.environ["BEGIN_TIME"],
    "end_time": os.environ["END_TIME"]
}, ensure_ascii=False))
'
    )"
    mcp_call_tool "chat_groups_list" "$args_json"
}

mcp_list_messages() {
    # Call chat_messages_list for a known chat_id. Useful for confirming that
    # a copied ID is valid before sending notifications to it.
    local chat_id="$1"
    local begin_time="${2:-$(default_begin_time)}"
    local end_time="${3:-$(default_end_time)}"
    local args_json

    require_cmd curl
    require_cmd python3

    args_json="$(
        CHAT_ID="$chat_id" BEGIN_TIME="$begin_time" END_TIME="$end_time" python3 -c '
import json
import os
print(json.dumps({
    "chat_id": os.environ["CHAT_ID"],
    "begin_time": os.environ["BEGIN_TIME"],
    "end_time": os.environ["END_TIME"]
}, ensure_ascii=False))
'
    )"
    mcp_call_tool "chat_messages_list" "$args_json"
}

send_wecom_message() {
    # metadata/tools.json says message_send requires chat_id, msg_type, and
    # text.content.
    # This function builds exactly that payload for plain text notifications.
    local message="${1:-job finished}"
    local message_json chat_id_json payload

    require_cmd curl
    require_cmd python3
    [[ -n "$wecom_chat_id" ]] || {
        echo "[ERROR] WECOM_CHAT_ID is required. Run 'bash scripts/wecom_mcp.sh list-tools' or check chat_groups_list/chat_messages_list output to get a chat_id." >&2
        return 1
    }

    message_json="$(printf '%s' "$message" | json_escape)"
    chat_id_json="$(printf '%s' "$wecom_chat_id" | json_escape)"

    payload="$(
        MESSAGE_JSON="$message_json" CHAT_ID_JSON="$chat_id_json" TOOL="$mcp_tool" python3 -c '
import json
import os

tool = os.environ["TOOL"]
chat_id = json.loads(os.environ["CHAT_ID_JSON"])
message = json.loads(os.environ["MESSAGE_JSON"])
payload = {
    "jsonrpc": "2.0",
    "id": 3,
    "method": "tools/call",
    "params": {
        "name": tool,
        "arguments": {
            "chat_id": chat_id,
            "msg_type": "text",
            "text": {
                "content": message
            }
        }
    }
}
print(json.dumps(payload, ensure_ascii=False))
'
    )"

    mcp_initialize
    mcp_post "$payload"
}

main() {
    # Convenience CLI: "send" is explicit, and any unknown first token is
    # treated as message text so `bash scripts/wecom_mcp.sh "hello"` still works.
    local command_name="${1:-send}"

    case "$command_name" in
        -h|--help|help)
            usage
            ;;
        list-tools)
            require_cmd curl
            require_cmd python3
            mcp_list_tools
            ;;
        list-groups)
            shift || true
            mcp_list_groups "${1:-}" "${2:-}"
            ;;
        list-messages)
            shift || true
            if [[ "${1:-}" == "--chat-id" ]]; then
                wecom_chat_id="$2"
                shift 2
            fi
            [[ -n "$wecom_chat_id" ]] || {
                echo "[ERROR] WECOM_CHAT_ID is required. Use --chat-id CHAT_ID or set WECOM_CHAT_ID." >&2
                return 1
            }
            mcp_list_messages "$wecom_chat_id" "${1:-}" "${2:-}"
            ;;
        send)
            shift || true
            if [[ "${1:-}" == "--chat-id" ]]; then
                wecom_chat_id="$2"
                shift 2
            fi
            send_wecom_message "${1:-job finished}"
            ;;
        *)
            send_wecom_message "$*"
            ;;
    esac
}

if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
    main "$@"
fi

"""Install via pip using IPv4 when this host's IPv6 downloads stall.

This affects only this installer process, never system network settings.
"""
import runpy
import socket
import sys

original = socket.getaddrinfo
def ipv4(host, port, family=0, type=0, proto=0, flags=0):
    return original(host, port, socket.AF_INET, type, proto, flags)
socket.getaddrinfo = ipv4
sys.argv[0] = "pip"
runpy.run_module("pip", run_name="__main__")

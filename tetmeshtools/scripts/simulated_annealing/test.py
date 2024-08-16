# import socket
# import ipaddress
# import re
# import subprocess
# # def send_wol_packet(mac_address, public_ip, port=9):
# #     # Convert MAC address to bytes
# #     mac_bytes = bytes.fromhex(mac_address.replace(':', '').replace('-', ''))
# #     # Create the magic packet
# #     magic_packet = b'\xff' * 6 + mac_bytes * 16
# #     # Send the packet to the public IP and port
# #     sock = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
# #     sock.sendto(magic_packet, (public_ip, port))

# # # Replace with your computer's MAC address and your public IP
# mac_address = '68-54-5A-CC-D6-1B'
# # public_ip = '10.68.142.181'

# def get_local_ip_address():
#     hostname = socket.gethostname()
#     local_ip = socket.gethostbyname(hostname)
#     return local_ip
# # send_wol_packet(mac_address, public_ip)

# def calculate_broadcast_address(ip, subnet_mask):
#     network = ipaddress.IPv4Network(f"{ip}/{subnet_mask}", strict=False)
#     return network.broadcast_address

# def send_wol_packet(mac_address, broadcast_ip='255.255.255.255', port=9):
#     mac_bytes = bytes.fromhex(mac_address.replace(':', '').replace('-', ''))
#     magic_packet = b'\xff' * 6 + mac_bytes * 16
#     sock = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
#     sock.setsockopt(socket.SOL_SOCKET, socket.SO_BROADCAST, 1)
#     sock.sendto(magic_packet, (broadcast_ip, port))

# # Example usage:
# # print(get_local_ip_address())
# ip_address = get_local_ip_address()
# broadcast_address = calculate_broadcast_address(ip_address, "255.255.240.0")
# broadcast_address = str(broadcast_address)
# # print(broadcast_address)
# broadcast_address = "10.68.138.129"
# send_wol_packet(mac_address, broadcast_ip=broadcast_address)

# def get_ip_from_mac(mac_address):
#     # Perform ARP command to fetch ARP table
#     arp_output = subprocess.check_output(['arp', '-a']).decode('utf-8')

#     # Search ARP table for the MAC address
#     pattern = r'([0-9]+\.[0-9]+\.[0-9]+\.[0-9]+)\s+([0-9a-fA-F]{2}-[0-9a-fA-F]{2}-[0-9a-fA-F]{2}-[0-9a-fA-F]{2}-[0-9a-fA-F]{2}-[0-9a-fA-F]{2})'
#     matches = re.findall(pattern, arp_output)
#     mac_address = mac_address.replace('-','').replace(':','').lower()
#     for match in matches:
#         full_mac = match[1].replace('-', '')
#         print(full_mac, mac_address)
#         if match[1].lower() == mac_address:
#             return match[0]

#     return None
    


# mac_address = "00-00-be-ef-ca-fe"
# print(get_ip_from_mac(mac_address))


# from wakeonlan import send_magic_packet
import socket
import struct

def mac_str_to_bytes(mac_str):
    # Remove any non-hexadecimal characters (like colons or hyphens)
    mac_hex = mac_str.replace(':', '').replace('-', '')
    # Convert the string to bytes
    return bytes.fromhex(mac_hex)

def wake_on_lan(mac_address, laptop_ip):
    # Convert MAC address string to bytes
    server_mac_address = mac_str_to_bytes(mac_address)
    
    # Create a raw socket
    with socket.socket(socket.AF_INET, socket.SOCK_RAW) as s:
        s.bind((laptop_ip, 0))  # Bind to any available interface

        # Set socket options to enable sending broadcast packet
        s.setsockopt(socket.SOL_SOCKET, socket.SO_BROADCAST, 1)

        # Construct the magic packet
        magic_packet = b'\xFF' * 6 + server_mac_address * 16

        # Send the magic packet
        s.sendto(magic_packet, ('<broadcast>', 7))  # Port number 7 is used for Wake-on-LAN


def get_ip_address():
    # Create a socket object
    s = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)

    try:
        # Connect to a remote server (doesn't have to be reachable)
        s.connect(('8.8.8.8', 80))  # Using Google's DNS server as an example
        # Get the local IP address bound to the socket
        ip_address = s.getsockname()[0]
    except Exception as e:
        print(f"Error: {e}")
        ip_address = None
    finally:
        s.close()

    return ip_address

# Example usage:
if __name__ == "__main__":
    server_mac = "68-54-5A-CC-D6-1B"  # Replace with your server's MAC address
    laptop_ip = get_ip_address()  # Replace with your laptop's IP address
    print(':'.join(format(x, '02x') for x in mac_str_to_bytes(server_mac)))
    wake_on_lan(server_mac, laptop_ip)


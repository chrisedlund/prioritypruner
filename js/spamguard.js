function spamGuard(name,domain,com) {
    var email = name + "@" + domain + "." + com;
    location = "mailto:" + email + "?subject=PriorityPruner";
}
